rm(list=ls())
#####
library(dismo)
library(rgeos)
library(mosaic)
library(dplyr)
library(geosphere)
library(maps)
library(raster)
library(sp)
library(gridExtra)
library(RColorBrewer)
library(dbscan)
library(beepr)
############ 
############ 
############ MUTASA ###############
############
############ 
##### read in data #####
setwd()
Mutasa4 <- read.csv("Mutasa4.csv", header=TRUE, sep=",")
##### remove those with less than 500 points #####
Mutasa4b <- Mutasa4[,-c(20:21)]

Mutasa4b <- droplevels(Mutasa4b[which(Mutasa4b$tot_points > 500),])
Mutasa4b <- Mutasa4b[order(Mutasa4b$partid, Mutasa4b$point_num),]

### fix those with excessive points taken affecting clusters (points every 10 seconds or less when 15+ in a row)--> make once a min #####
parts_too_close <- as.factor(names(which(table(as.factor(Mutasa4b$partid[which((Mutasa4b$time_atpoint) < 0.00011574)])) > 10)))
for(ii in 1:length(parts_too_close)){
  ind_keep <- NA
  ind_rem <- NA
  part <- parts_too_close[ii]
  temp <- Mutasa4b[which(Mutasa4b$partid == part),]
  inds <- which(diff(temp$log_date_time_num) < 0.00011574)
  counts <- cumsum(c(1,diff(inds)!=1))
  counts[which(counts %in% c(which(summary(as.factor(counts), maxsum=1000) < 15)))]<- NA
  counts <- as.factor(counts)
  if(nlevels(counts) > 0){
    for(jj in 1:nlevels(counts)){
      count <- levels(counts)[jj]
      inds2 <- inds[which(counts == count)]
      ind_keep <- c(ind_keep,inds2[seq(from=1,to=length(inds2), by=12)])
      ind_rem <- c(ind_rem,inds2[-c(seq(from=1,to=length(inds2), by=12))])
    }
    ind_rem <- ind_rem[-1]
    Mutasa4b <- Mutasa4b[-c(which(Mutasa4b$partid == part)[ind_rem]),]
    temp <- Mutasa4b[which(Mutasa4b$partid == part),]
    time_at_points <- c(c(diff(c(temp$log_date_time_num[1], temp$log_date_time_num[2]))/2),
                        c(diff(c(temp$log_date_time_num),2)/2),
                        c(diff(c(temp$log_date_time_num[c(nrow(temp)-1)], temp$log_date_time_num[c(nrow(temp))]))/2))
    Mutasa4b$time_atpoint[which(Mutasa4b$partid == part)] <- time_at_points
    Mutasa4b$total_time[which(Mutasa4b$partid == part)] <- sum(time_at_points)
  }
}
##
#### define clusters with HDBSCAN and DBSCAN #####
Mutasa4b$clust <- numeric(length=nrow(Mutasa4b))
Mutasa4b$points_inclust <- numeric(length=nrow(Mutasa4b))
Mutasa4b$time_inclust <- numeric(length=nrow(Mutasa4b))

parts <- levels(as.factor(Mutasa4b$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Mutasa4c <- Mutasa4b[which(Mutasa4b$partid == part),]
  xy <- SpatialPointsDataFrame(
    matrix(c(Mutasa4c$log_long, Mutasa4c$log_lat), ncol=2), data.frame(ID=seq(1:nrow(Mutasa4c))),
    proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  mdist <- distm(xy)
  bb <- hdbscan(as.dist(mdist), minPts = 15) 
  ##  if more than 2 clusters (and noise), all good
  if(length(unique(bb$cluster)) > 3){
    xy$clust <- bb$cluster
  }
  ## if 2 clusters (and noise)
  else if(length(unique(bb$cluster)) == 3){
    res <- .Call('_dbscan_computeStability', bb$hc, 15, compute_glosh=TRUE)
    if(res$'0'$stability < (res$'1'$stability + res$'2'$stability)){ ## 2 clusters are good, keep going as normal
      xy$clust <- bb$cluster
    }
    else{ ## root is better --> still want to find outliers though
      res2 <- dbscan(as.dist(mdist), eps = 25, minPts = 6)
      clust_temp <- as.numeric(as.character(as.factor(res2$cluster)[(which.max(summary(as.factor(res2$cluster))))])) ## not actual cluster # with which.max since minPts > 1
      xy$clust <- 0
      xy$clust[which(res2$cluster == clust_temp)] <- 1
    }
  }
  else if(length(unique(bb$cluster)) == 2){ ## either 1 clust and outliers or 2 clusts and no outliers; either way, should be good
    xy$clust <- bb$cluster
  }
  else if(length(unique(bb$cluster)) == 1){ ## only 1 clust found but all considered 'outliers'
    res2 <- dbscan(as.dist(mdist), eps = 25, minPts = 1)
    clust_temp <- which.max(summary(as.factor(res2$cluster)))
    xy$clust <- 0
    xy$clust[which(res2$cluster == clust_temp)] <- 1
  }
  if(length(which(bb$cluster == 0))/length(bb$cluster) > 0.5){
    res2 <- dbscan(as.dist(mdist), eps = 25, minPts = 1)
    clust_temp <- which.max(summary(as.factor(res2$cluster))) ## actual clust # since no noise points (since minPts =1)
    xy$clust <- 0
    xy$clust[which(res2$cluster == clust_temp)] <- 1
  }
  xy2 <- as.data.frame(xy)
  Mutasa4c$clust <- factor(as.numeric(xy2$clust))
  xx <- aggregate(time_atpoint ~ clust, Mutasa4c, sum)
  for(jj in 1:nrow(xx)){
    clust <- xx$clust[jj]
    
    Mutasa4c$points_inclust[which(Mutasa4c$clust == clust)] <- length(which(xy2$clust == clust))
    Mutasa4c$time_inclust[which(Mutasa4c$clust == clust)] <- xx$time_atpoint[jj]
  }
  Mutasa4b$clust[which(Mutasa4b$partid == part)] <- as.numeric(as.character(Mutasa4c$clust))
  Mutasa4b$points_inclust[which(Mutasa4b$partid == part)] <- Mutasa4c$points_inclust
  Mutasa4b$time_inclust[which(Mutasa4b$partid == part)] <- Mutasa4c$time_inclust
}
Mutasa4b$time_inclust_hour <- Mutasa4b$time_inclust*24
Mutasa4b$time_inclust_min <- Mutasa4b$time_inclust*24*60
Mutasa_short1 <- Mutasa4b[!duplicated(Mutasa4b[,c(1,21)]),]

Mutasa_short <- Mutasa_short1[which(Mutasa_short1$clust != 0),]
#beep()

### remove points greater than 1000 meters from center of cluster ####
Mutasa4b_short <- Mutasa4b[which(Mutasa4b$clust != 0),]
parts <- levels(as.factor(Mutasa4b_short$partid))
for(ii in 1:nlevels(as.factor(Mutasa4b_short$partid))){
  part <- parts[ii]
  temp <- Mutasa4b_short[which(Mutasa4b_short$partid == part),]
  out <- tapply(1:nrow(temp), as.factor(temp$clust), function(x) max(distm(as.matrix(temp[x,10:9], fun=distHaversine)))) ## meters
  clust_cut <- which(out > 2000)
  if(length(clust_cut) > 0){
    cluster_centers <- aggregate(temp[,10:9], list(Clust = as.factor(temp$clust)), median)
    for(jj in 1:length(clust_cut)){
      clust <- clust_cut[jj]
      temp2 <- temp[which(temp$clust == clust),]
      out2 <- distm(as.matrix(temp2[,10:9]), cluster_centers[which(cluster_centers$Clust==clust),2:3], fun=distHaversine)
      cuts <- which(out2 > 1000)
      Mutasa4b_short$clust[c(c(which(Mutasa4b_short$partid == part & Mutasa4b_short$clust == clust))[cuts])] <- 0
      Mutasa4b$clust[c(c(which(Mutasa4b$partid == part & Mutasa4b$clust == clust))[cuts])] <- 0
      
      Mutasa4b_short$points_inclust[c(which(Mutasa4b_short$partid == part & Mutasa4b_short$clust == clust))] <- nrow(Mutasa4b_short[c(which(Mutasa4b_short$partid == part & Mutasa4b_short$clust == clust)),])
      Mutasa4b$points_inclust[c(which(Mutasa4b$partid == part & Mutasa4b$clust == clust))] <- nrow(Mutasa4b[c(which(Mutasa4b$partid == part & Mutasa4b$clust == clust)),])
    }
  }
}
Mutasa4b_short2 <- Mutasa4b_short[which(Mutasa4b_short$clust != 0),]
Mutasa_short <- Mutasa4b_short2[!duplicated(Mutasa4b_short2[,c(1,21)]),]
##
###### limit number of clusters based on time spent and points ####
## by time spent
Mutasa_short_time <- Mutasa_short[order(Mutasa_short$partid, -Mutasa_short$time_inclust),]
### cut at 20 locations per person
Mutasa_short_time$in_out <- factor(NA, levels=c(0,1))
for(ii in 1:length(parts)){
  part <- parts[ii]
  ind <- which(Mutasa_short_time$partid == part)
  Mutasa_short_time$in_out[ind[1:20]] <- 1
  if(length(ind) > 20){
    Mutasa_short_time$in_out[ind[21:length(ind)]] <- 0
  }
}

## by number of points
Mutasa_short_points <- Mutasa_short[order(Mutasa_short$partid, -Mutasa_short$points_inclust),]
Mutasa_short_points$in_out <- factor(NA, levels=c(0,1))
for(ii in 1:length(parts)){
  part <- parts[ii]
  ind <- which(Mutasa_short_points$partid == part)
  Mutasa_short_points$in_out[ind[1:20]] <- 1
  if(length(ind) > 20){
    Mutasa_short_points$in_out[ind[21:length(ind)]] <- 0
  }
}
##
parts <- levels(as.factor(Mutasa_short$partid))
counts <- numeric(length=length(parts))
for(ii in 1:length(parts)){
  part <- parts[ii]
  counts[ii] <- length(sort(unique(c(Mutasa_short_points$clust[which(Mutasa_short_points$partid == part & Mutasa_short_points$in_out == 1)],
                                     Mutasa_short_time$clust[which(Mutasa_short_time$partid == part & Mutasa_short_time$in_out == 1)]))))
}

Mutasa_short_points2 <- Mutasa_short_points[which(Mutasa_short_points$in_out == 1),]
Mutasa_short_time2 <- Mutasa_short_time[which(Mutasa_short_time$in_out == 1),]

##### remove clusters with less than 30 minutes or less than 5 points #####
Mutasa_short2 <- Mutasa_short
Mutasa_short2$keep <- factor(0, levels=c(0,1))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- sort(unique(c(Mutasa_short_points$clust[which(Mutasa_short_points$partid == part & Mutasa_short_points$in_out == 1)],
                          Mutasa_short_time$clust[which(Mutasa_short_time$partid == part & Mutasa_short_time$in_out == 1)])))
  Mutasa_short2$keep[which(Mutasa_short2$partid == part & Mutasa_short2$clust %in% clusts)] <- 1
}

Mutasa_short2b <- Mutasa_short2[which(Mutasa_short2$keep == 1),]

Mutasa_short2c <- Mutasa_short2b
a <- length(unique(Mutasa_short2c$partid[which(Mutasa_short2c$time_inclust_min < 30 | Mutasa_short2c$points_inclust < 5)])) ## 1
if(a == 0){
  Mutasa_short2d <- Mutasa_short2c}
if(a > 0){
  Mutasa_short2d <- Mutasa_short2c[-c(which(Mutasa_short2c$time_inclust_min < 30 | Mutasa_short2c$points_inclust < 5)),]
}

#####
Mutasa5 <- Mutasa4b
Mutasa5$keep <- factor(0, levels=c(0,1))
parts <- levels(as.factor(Mutasa5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- Mutasa_short2d$clust[which(Mutasa_short2d$partid == part)]
  Mutasa5$keep[which(Mutasa5$partid == part & Mutasa5$clust %in% clusts)] <- 1
}

### define home cluster by distance from GPS of home listed #####
Mutasa_short2d$home_clust <- factor(0, levels=c(0,1))
Mutasa5$home_clust <- factor(0, levels=c(0,1))

parts <- levels(as.factor(Mutasa5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  min_dist <- min(Mutasa_short2d$hhdist_m_haversine[which(Mutasa_short2d$partid == part)])
  home_clust <- Mutasa_short2d$clust[which(Mutasa_short2d$partid == part & Mutasa_short2d$hhdist_m_haversine == min_dist)]
  Mutasa_short2d$home_clust[which(Mutasa_short2d$partid == part & Mutasa_short2d$clust == home_clust)] <- 1
  Mutasa5$home_clust[which(Mutasa5$partid == part & Mutasa5$clust == home_clust)] <- 1
}
#
##### define trips to locations #####
Mutasa5$time_atpoint_hour <- Mutasa5$time_atpoint*24
Mutasa5$time_atpoint_min <- Mutasa5$time_atpoint*24*60

Mutasa5$trip_count_new <- numeric(length=nrow(Mutasa5))
Mutasa5$trip_time_new <- numeric(length=nrow(Mutasa5))
parts <- levels(as.factor(Mutasa5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Mutasa5$clust[which(Mutasa5$partid == part & Mutasa5$keep == 1)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    inds <- which(Mutasa5$partid == part & Mutasa5$clust == clust)
    counts <- cumsum(c(1,diff(inds)>6))
    
    counts[which(counts %in% c(which(summary(as.factor(counts), maxsum=500) < 7)))]<- NA
    counts <- as.factor(counts)
    levels(counts) <- c(1:nlevels(counts))
    Mutasa5$trip_count_new[which(Mutasa5$partid == part & Mutasa5$clust == clust)] <- as.numeric(as.character(counts))
    Mutasa5c <- Mutasa5[which(Mutasa5$partid == part & Mutasa5$clust == clust),]
    if(length(which(!(is.na(Mutasa5c$trip_count_new))))){
      xx <- aggregate(time_atpoint ~ trip_count_new, Mutasa5c, sum)
      for(kk in 1:nrow(xx)){
        Mutasa5$trip_time_new[which(Mutasa5$partid == part & Mutasa5$clust == clust & Mutasa5$trip_count_new == kk)] <- xx$time_atpoint[kk]
      }
    }
  }
}
Mutasa5$trip_time_new[which(is.na(Mutasa5$trip_count_new))] <- NA
Mutasa5$trip_time_new_hour <- Mutasa5$trip_time_new*24
Mutasa5$trip_time_new_min <- Mutasa5$trip_time_new*24*60
#beep()
##### redfeine home cluster ###### 
### do 5 closest distances --> only if less than 1-2 km --> which spend most time at --> make that home cluster
Mutasa5$home_clust2 <- factor(0, levels=c(0,1))

parts <- levels(as.factor(Mutasa5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Mutasa5aa <- Mutasa5[which(Mutasa5$partid==part),]
  Mutasa5aa <- Mutasa5aa[order(Mutasa5aa$hhdist_m_haversine),]
  Mutasa5aa2 <- Mutasa5aa[!duplicated(Mutasa5aa$clust),]
  min_dists <- sort(Mutasa5aa2$hhdist_m_haversine[which(Mutasa5aa2$clust != 0 & Mutasa5aa2$keep==1)])
  min_dists2 <- min_dists[which(min_dists < 500)]
  clust_temp <- Mutasa5[which(Mutasa5$partid == part & Mutasa5$hhdist_m_haversine %in% min_dists2),]
  home_clust2 <- clust_temp$clust[which.max(clust_temp$time_inclust)]
  Mutasa5$home_clust2[which(Mutasa5$partid == part & Mutasa5$clust == home_clust2)] <- 1
}
#
### add one point per hour in location if <9 hrs missing overnight or <4hrs missing during the same day and points before and after gap in same cluster ####
Mutasa5$date_time_new <- as.POSIXct((Mutasa5$log_date_time_num*24*60*60), tz="Africa/Lusaka", origin="1899-12-30")
dates <- c(unlist(strsplit(as.character(Mutasa5$date_time_new)," "))[c(seq(from=1,to=c((nrow(Mutasa5)*2)-1), by=2))])
times <- c(unlist(strsplit(as.character(Mutasa5$date_time_new)," "))[c(seq(from=2,to=c(nrow(Mutasa5)*2), by=2))])
hours <- c(unlist(strsplit(times,":"))[c(seq(from=1,to=c((nrow(Mutasa5)*3)-2), by=3))])
Mutasa5$date_new <- dates
Mutasa5$time_new <- times
Mutasa5$time_new_hour <- as.numeric(as.character(hours))

Mutasa5_temp <- Mutasa5[order(Mutasa5$partid, Mutasa5$date_time_new),]

times <- 0
times2 <- 0
times_all <- 0

parts <- levels(as.factor(Mutasa5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Mutasa5_temp$clust[which(Mutasa5_temp$partid == part & Mutasa5_temp$keep == 1)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    temp <- Mutasa5_temp[which(Mutasa5_temp$partid == part & Mutasa5_temp$clust == clust & Mutasa5_temp$keep == 1),c(36,38)]
    for(kk in 1:(nrow(temp)-1)){
      if((as.Date(temp$date_new[kk])+1) == temp$date_new[kk + 1]){
        if(temp$time_new_hour[kk] != 23 | temp$time_new_hour[kk + 1] != 0){
          if((as.numeric(rownames(temp))[kk] + 1) == as.numeric(rownames(temp))[kk + 1]){
            length(times) <- 0
            length(times2) <- 0
            length(times_all) <- 0
            if((temp$time_new_hour[kk] != 23) & (temp$time_new_hour[kk + 1] != 0)){
              times <- seq((temp$time_new_hour[kk]+1),23,by=1)
              times2 <- seq(0,(temp$time_new_hour[kk + 1] - 1),by=1)
              times_all <- c(times,times2)
              add <- Mutasa5_temp[which(Mutasa5_temp$partid == part & Mutasa5_temp$clust == clust  & Mutasa5_temp$keep == 1)[kk],]
              add[,c(4:8,11:13,22:25,28:29,31:33,35)] <- NA
              adds <- as.data.frame(lapply(add, rep, length(times_all)))
              adds$date_new <- c(rep(temp$date_new[kk], times=length(times)), rep(temp$date_new[kk + 1], times=length(times2)))
              adds$time_new_hour <- times_all
              adds$time_new <- paste0(times_all, ":00:00")
            }
            if((temp$time_new_hour[kk] != 23) & (temp$time_new_hour[kk + 1] == 0)){
              times <- seq((temp$time_new_hour[kk]+1),23,by=1)
              times_all <- times
              add <- Mutasa5_temp[which(Mutasa5_temp$partid == part & Mutasa5_temp$clust == clust  & Mutasa5_temp$keep == 1)[kk],]
              add[,c(4:8,11:13,22:25,28:29,31:33,35)] <- NA
              adds <- as.data.frame(lapply(add, rep, length(times_all)))
              adds$date_new <- c(rep(temp$date_new[kk], times=length(times)))
              adds$time_new_hour <- times_all
              adds$time_new <- paste0(times_all, ":00:00")
            }
            if((temp$time_new_hour[kk] == 23) & (temp$time_new_hour[kk + 1] != 0)){
              times2 <- seq(0,(temp$time_new_hour[kk + 1] - 1),by=1)
              times_all <- times2
              add <- Mutasa5_temp[which(Mutasa5_temp$partid == part & Mutasa5_temp$clust == clust  & Mutasa5_temp$keep == 1)[kk],]
              add[,c(4:8,11:13,22:25,28:29,31:33,35)] <- NA
              adds <- as.data.frame(lapply(add, rep, length(times_all)))
              adds$date_new <- c(rep(temp$date_new[kk + 1], times=length(times2)))
              adds$time_new_hour <- times_all
              adds$time_new <- paste0(times_all, ":00:00")
            }
            if(length(times_all) < 9){
              Mutasa5_temp <- rbind(Mutasa5_temp,adds)
            }
          }
        }
      }
      if((temp$date_new[kk] == temp$date_new[kk + 1]) & (temp$time_new_hour[kk] != temp$time_new_hour[kk + 1])){
        if((temp$time_new_hour[kk] + 1) != temp$time_new_hour[kk + 1]){
          if((as.numeric(rownames(temp))[kk] + 1) == as.numeric(rownames(temp))[kk + 1]){
            add <- Mutasa5_temp[which(Mutasa5_temp$partid == part & Mutasa5_temp$clust == clust  & Mutasa5_temp$keep == 1)[kk],]
            add[,c(4:8,11:13,22:25,28:29,31:33,35)] <- NA
            count <- (temp$time_new_hour[kk + 1] - temp$time_new_hour[kk])-1
            if(count < 4){
              count2 <- seq((temp$time_new_hour[kk] + 1), by=1, length.out=count)
              adds <- as.data.frame(lapply(add, rep, count))
              adds$time_new_hour <- count2
              adds$time_new <- paste0(count2, ":00:00")
              Mutasa5_temp <- rbind(Mutasa5_temp,adds)
            }
          }
        }
      }
    }
  }
}

Mutasa5_temp$time_of_day <- character(length=nrow(Mutasa5_temp))
Mutasa5_temp$time_of_day[which(Mutasa5_temp$time_new_hour >= 6 & Mutasa5_temp$time_new_hour <= 11)] <- "morning"  ## 6am-11am
Mutasa5_temp$time_of_day[which(Mutasa5_temp$time_new_hour >= 12 & Mutasa5_temp$time_new_hour <= 17)] <- "midday"  ## 12pm-5pm
Mutasa5_temp$time_of_day[which(Mutasa5_temp$time_new_hour >= 18 & Mutasa5_temp$time_new_hour <= 23)] <- "evening"  ## 6pm-11pm
Mutasa5_temp$time_of_day[which(Mutasa5_temp$time_new_hour <= 5)] <- "night"  ## 12am-5am
Mutasa5_temp$time_of_day <- factor(Mutasa5_temp$time_of_day, levels=c("morning","midday","evening","night"))

Mutasa5b <- droplevels(Mutasa5_temp[which(Mutasa5_temp$keep == 1 & !(is.na(Mutasa5_temp$trip_count_new))),])

##
###### calculate new values for time per trip and per cluster #####
Mutasa5b2 <- Mutasa5b[order(Mutasa5b$partid, Mutasa5b$date_new, Mutasa5b$time_new_hour, Mutasa5b$time_new),]

Mutasa5b2$time_atpoint_new <- as.numeric(NA)
Mutasa5b2$trip_time_new2 <- as.numeric(NA)
Mutasa5b2$trip_time_new2_morning <- as.numeric(NA)
Mutasa5b2$trip_time_new2_midday <- as.numeric(NA)
Mutasa5b2$trip_time_new2_evening <- as.numeric(NA)
Mutasa5b2$trip_time_new2_night <- as.numeric(NA)
Mutasa5b2$clust_trips_time_new <- as.numeric(NA)
Mutasa5b2$clust_trips_time_new_morning <- as.numeric(NA)
Mutasa5b2$clust_trips_time_new_midday <- as.numeric(NA)
Mutasa5b2$clust_trips_time_new_evening <- as.numeric(NA)
Mutasa5b2$clust_trips_time_new_night <- as.numeric(NA)

parts <- levels(as.factor(Mutasa5b2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Mutasa5b2$clust[which(Mutasa5b2$partid == part)]))
  for(jj in 1:length(clusts)){
    clust <- clusts[jj]
    trips <- levels(as.factor(Mutasa5b2$trip_count_new[which(Mutasa5b2$partid == part & Mutasa5b2$clust == clust)]))
    for(kk in 1:length(trips)){
      trip <- trips[kk]
      temp2 <- Mutasa5b2[which(Mutasa5b2$partid == part & Mutasa5b2$clust == clust & Mutasa5b2$trip_count_new == trip),]
      time_at_points_middle <- diff(c(as.POSIXct(paste(temp2$date_new, temp2$time_new), format="%Y-%m-%d %H:%M:%S", tz="Africa/Lusaka")),2)
      if(units(time_at_points_middle) == "hours"){
        time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/24
      }
      if(units(time_at_points_middle) == "mins"){
        time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/60/24
      }
      if(units(time_at_points_middle) == "secs"){
        time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/60/60/24
      }
      time_at_points <- c(as.numeric(difftime(as.POSIXct(paste(temp2$date_new[2], temp2$time_new[2])), as.POSIXct(paste(temp2$date_new[1], temp2$time_new[1])), unit="min"))/2/60/24, 
                          time_at_points_middle2,
                          as.numeric(difftime(as.POSIXct(paste(temp2$date_new[c(nrow(temp2))], temp2$time_new[c(nrow(temp2))])), 
                                              as.POSIXct(paste(temp2$date_new[c(nrow(temp2)-1)], temp2$time_new[c(nrow(temp2)-1)])), 
                                              unit="min"))/2/60/24)
      Mutasa5b2$time_atpoint_new[which(Mutasa5b2$partid == part & Mutasa5b2$clust == clust & Mutasa5b2$trip_count_new == trip)] <- time_at_points
    }
    temp3 <- Mutasa5b2[which(Mutasa5b2$partid == part & Mutasa5b2$clust == clust),]
    xx <- aggregate(time_atpoint_new ~ trip_count_new, temp3, sum)
    for(kk in 1:nrow(xx)){
      Mutasa5b2$trip_time_new2[which(Mutasa5b2$partid == part & Mutasa5b2$clust == clust & Mutasa5b2$trip_count_new == kk)] <- xx$time_atpoint[kk]
    }
    xx2 <- aggregate(time_atpoint_new ~ trip_count_new + time_of_day, temp3, sum)
    xx2_morning <- xx2[which(xx2$time_of_day == "morning"),]
    xx2_midday <- xx2[which(xx2$time_of_day == "midday"),]
    xx2_evening <- xx2[which(xx2$time_of_day == "evening"),]
    xx2_night <- xx2[which(xx2$time_of_day == "night"),]
    for(kk in 1:nrow(xx2_morning)){
      Mutasa5b2$trip_time_new2_morning[which(Mutasa5b2$partid == part & Mutasa5b2$clust == clust & Mutasa5b2$trip_count_new == xx2_morning$trip_count_new[kk])] <- xx2_morning$time_atpoint_new[kk]
    }
    for(kk in 1:nrow(xx2_midday)){
      Mutasa5b2$trip_time_new2_midday[which(Mutasa5b2$partid == part & Mutasa5b2$clust == clust & Mutasa5b2$trip_count_new == xx2_midday$trip_count_new[kk])] <- xx2_midday$time_atpoint_new[kk]
    }
    for(kk in 1:nrow(xx2_evening)){
      Mutasa5b2$trip_time_new2_evening[which(Mutasa5b2$partid == part & Mutasa5b2$clust == clust & Mutasa5b2$trip_count_new == xx2_evening$trip_count_new[kk])] <- xx2_evening$time_atpoint_new[kk]
    }
    for(kk in 1:nrow(xx2_night)){
      Mutasa5b2$trip_time_new2_night[which(Mutasa5b2$partid == part & Mutasa5b2$clust == clust & Mutasa5b2$trip_count_new == xx2_night$trip_count_new[kk])] <- xx2_night$time_atpoint_new[kk]
    }
  }
  temp4 <- Mutasa5b2[which(Mutasa5b2$partid == part),]
  temp4b <- temp4[!duplicated(temp4[c(1,21,30)]),]
  xx <- aggregate(trip_time_new2 ~ clust, temp4b, sum)
  for(kk in 1:nrow(xx)){
    Mutasa5b2$clust_trips_time_new[which(Mutasa5b2$partid == part & Mutasa5b2$clust == xx$clust[kk] & !(is.na(Mutasa5b2$trip_count_new)))] <- xx$trip_time_new2[kk]
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_morning))),])>0){
    morning <- aggregate(trip_time_new2_morning ~ clust, temp4b, sum)
    if(nrow(morning) > 0){
      for(kk in 1:nrow(morning)){
        Mutasa5b2$clust_trips_time_new_morning[which(Mutasa5b2$partid == part & Mutasa5b2$clust == morning$clust[kk] & !(is.na(Mutasa5b2$trip_count_new)))] <- morning$trip_time_new2_morning[kk]
      }
    }
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_midday))),])>0){
    midday <- aggregate(trip_time_new2_midday ~ clust, temp4b, sum)
    if(nrow(midday) > 0){
      for(kk in 1:nrow(midday)){
        Mutasa5b2$clust_trips_time_new_midday[which(Mutasa5b2$partid == part & Mutasa5b2$clust == midday$clust[kk] & !(is.na(Mutasa5b2$trip_count_new)))] <- midday$trip_time_new2_midday[kk]
      }
    }
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_evening))),])>0){
    evening <- aggregate(trip_time_new2_evening ~ clust, temp4b, sum)
    if(nrow(evening) > 0){
      for(kk in 1:nrow(evening)){
        Mutasa5b2$clust_trips_time_new_evening[which(Mutasa5b2$partid == part & Mutasa5b2$clust == evening$clust[kk] & !(is.na(Mutasa5b2$trip_count_new)))] <- evening$trip_time_new2_evening[kk]
      }
    }
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_night))),])>0){
    night <- aggregate(trip_time_new2_night ~ clust, temp4b, sum)
    if(nrow(night) > 0){
      for(kk in 1:nrow(night)){
        Mutasa5b2$clust_trips_time_new_night[which(Mutasa5b2$partid == part & Mutasa5b2$clust == night$clust[kk] & !(is.na(Mutasa5b2$trip_count_new)))] <- night$trip_time_new2_night[kk]
      }
    }
  }
}


Mutasa5b2a <- Mutasa5b2[order(Mutasa5b2$partid, Mutasa5b2$clust, Mutasa5b2$trip_count_new),]
Mutasa5b2a2 <- Mutasa5b2a[which(!(is.na(Mutasa5b2a$trip_time_new2))),]


Mutasa5b2a2$part_trips_time_new <- as.numeric(0)
Mutasa5b2a2$part_trips_time_new_morning <- as.numeric(0)
Mutasa5b2a2$part_trips_time_new_midday <- as.numeric(0)
Mutasa5b2a2$part_trips_time_new_evening <- as.numeric(0)
Mutasa5b2a2$part_trips_time_new_night <- as.numeric(0)
temp <- Mutasa5b2a2[!duplicated(Mutasa5b2a2[c(1,21)]),]
tot <- aggregate(clust_trips_time_new ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Mutasa5b2a2$part_trips_time_new[which(Mutasa5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new[ii]
}
tot <- aggregate(clust_trips_time_new_morning ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Mutasa5b2a2$part_trips_time_new_morning[which(Mutasa5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_morning[ii]
}
tot <- aggregate(clust_trips_time_new_midday ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Mutasa5b2a2$part_trips_time_new_midday[which(Mutasa5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_midday[ii]
}
tot <- aggregate(clust_trips_time_new_evening ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Mutasa5b2a2$part_trips_time_new_evening[which(Mutasa5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_evening[ii]
}
tot <- aggregate(clust_trips_time_new_night ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Mutasa5b2a2$part_trips_time_new_night[which(Mutasa5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_night[ii]
}

#######
Mutasa5_new <- Mutasa5b2a2[,c(1:3,14,15,4,12,5:11,29,16:20,27,21:23,25,30,31,33,34,35:55,26)]

### take out those with <15 minutes in cluster or <15 minute trips that may be here now because of clusters that were >1000 meter radius ####
Mutasa5_new2 <- Mutasa5_new[which(Mutasa5_new$clust_trips_time_new > 0.01),]
Mutasa5_new2b <- Mutasa5_new2[which(Mutasa5_new2$trip_time_new2 > 0.01),]

Mutasa5_clusts_shorter <- Mutasa5_new2b[order(Mutasa5_new2b$partid, Mutasa5_new2b$clust),]
Mutasa5_clusts_shorter2 <- Mutasa5_clusts_shorter[!duplicated(Mutasa5_clusts_shorter[,c(1,22)]),c(1:7,16,17,12,13,19:25,29,30:35,46:50,41:45,51)]

Mutasa5_clust_trips_shorter <- Mutasa5_new2b[order(Mutasa5_new2b$partid, Mutasa5_new2b$clust, Mutasa5_new2b$trip_count),]
Mutasa5_clust_trips_shorter2 <- Mutasa5_clust_trips_shorter[!duplicated(Mutasa5_clust_trips_shorter[,c(1,22,26)]),c(1:7,16,17,12,13,19:22,24,26:28,29,30:35,46:50,41:45,36:40,51)]
##
#####
setwd()
write.csv(Mutasa4, "Mutasa4.csv", row.names = FALSE)
write.csv(Mutasa4b, "Mutasa4b.csv", row.names = FALSE)
write.csv(Mutasa5, "Mutasa5.csv", row.names = FALSE)
write.csv(Mutasa5_new2b, "Mutasa5_new.csv", row.names = FALSE)

write.csv(Mutasa5_clusts_shorter2, "Mutasa5_clusts_shorter2.csv", row.names = FALSE)
write.csv(Mutasa5_clust_trips_shorter2, "Mutasa5_clust_trips_shorter2.csv", row.names = FALSE)
##
############ 
############ 
############ CHOMA ###############
############ 
############ 
##### read in data ####
setwd()
Choma4 <- read.csv("Choma4.csv", header=TRUE, sep=",")
##### remove those with less than 500 points #####
Choma4b <- Choma4[,-c(20:21)]

Choma4b <- droplevels(Choma4b[which(Choma4b$tot_points > 500),])
Choma4b <- Choma4b[order(Choma4b$partid, Choma4b$point_num),]


### fix those with excessive points taken affecting clusters (points every 10 seconds or less when 15+ in a row)--> make once a min #####
parts_too_close <- as.factor(names(which(table(as.factor(Choma4b$partid[which((Choma4b$time_atpoint) < 0.00011574)])) > 10)))
for(ii in 1:length(parts_too_close)){
  ind_keep <- NA
  ind_rem <- NA
  part <- parts_too_close[ii]
  temp <- Choma4b[which(Choma4b$partid == part),]
  inds <- which(diff(temp$log_date_time_num) < 0.00011574)
  counts <- cumsum(c(1,diff(inds)!=1))
  counts[which(counts %in% c(which(summary(as.factor(counts), maxsum=1000) < 15)))]<- NA
  counts <- as.factor(counts)
  for(jj in 1:nlevels(counts)){
    count <- levels(counts)[jj]
    inds2 <- inds[which(counts == count)]
    ind_keep <- c(ind_keep,inds2[seq(from=1,to=length(inds2), by=12)])
    ind_rem <- c(ind_rem,inds2[-c(seq(from=1,to=length(inds2), by=12))])
  }
  ind_rem <- ind_rem[-1]
  Choma4b <- Choma4b[-c(which(Choma4b$partid == part)[ind_rem]),]
  temp <- Choma4b[which(Choma4b$partid == part),]
  time_at_points <- c(c(diff(c(temp$log_date_time_num[1], temp$log_date_time_num[2]))/2),
                      c(diff(c(temp$log_date_time_num),2)/2),
                      c(diff(c(temp$log_date_time_num[c(nrow(temp)-1)], temp$log_date_time_num[c(nrow(temp))]))/2))
  Choma4b$time_atpoint[which(Choma4b$partid == part)] <- time_at_points
  Choma4b$total_time[which(Choma4b$partid == part)] <- sum(time_at_points)
}

#### define clusters with HDBSCAN and DBSCAN #####
Choma4b$clust <- numeric(length=nrow(Choma4b))
Choma4b$points_inclust <- numeric(length=nrow(Choma4b))
Choma4b$time_inclust <- numeric(length=nrow(Choma4b))

parts <- levels(as.factor(Choma4b$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Choma4c <- Choma4b[which(Choma4b$partid == part),]
  xy <- SpatialPointsDataFrame(
    matrix(c(Choma4c$log_long, Choma4c$log_lat), ncol=2), data.frame(ID=seq(1:nrow(Choma4c))),
    proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  mdist <- distm(xy)
  bb <- hdbscan(as.dist(mdist), minPts = 15) 
  ##  if more than 2 clusters (and noise), all good
  if(length(unique(bb$cluster)) > 3){
    xy$clust <- bb$cluster
  }
  ## if 2 clusters (and noise)
  else if(length(unique(bb$cluster)) == 3){
    res <- .Call('_dbscan_computeStability', bb$hc, 15, compute_glosh=TRUE)
    if(res$'0'$stability < (res$'1'$stability + res$'2'$stability)){ ## 2 clusters are good, keep going as normal
      xy$clust <- bb$cluster
    }
    else{ ## root is better --> still want to find outliers though
      res2 <- dbscan(as.dist(mdist), eps = 25, minPts = 6)
      clust_temp <- as.numeric(as.character(as.factor(res2$cluster)[(which.max(summary(as.factor(res2$cluster))))])) ## not actual cluster # with which.max since minPts > 1
      xy$clust <- 0
      xy$clust[which(res2$cluster == clust_temp)] <- 1
    }
  }
  else if(length(unique(bb$cluster)) == 2){ ## either 1 clust and outliers or 2 clusts and no outliers; either way, should be good
    xy$clust <- bb$cluster
  }
  else if(length(unique(bb$cluster)) == 1){ ## only 1 clust found but all considered 'outliers'
    res2 <- dbscan(as.dist(mdist), eps = 25, minPts = 1)
    clust_temp <- which.max(summary(as.factor(res2$cluster)))
    xy$clust <- 0
    xy$clust[which(res2$cluster == clust_temp)] <- 1
  }
  if(length(which(bb$cluster == 0))/length(bb$cluster) > 0.5){
    res2 <- dbscan(as.dist(mdist), eps = 25, minPts = 1)
    clust_temp <- which.max(summary(as.factor(res2$cluster))) ## actual clust # since no noise points (since minPts =1)
    xy$clust <- 0
    xy$clust[which(res2$cluster == clust_temp)] <- 1
  }
  xy2 <- as.data.frame(xy)
  Choma4c$clust <- factor(as.numeric(xy2$clust))
  xx <- aggregate(time_atpoint ~ clust, Choma4c, sum)
  for(jj in 1:nrow(xx)){
    clust <- xx$clust[jj]
    
    Choma4c$points_inclust[which(Choma4c$clust == clust)] <- length(which(xy2$clust == clust))
    Choma4c$time_inclust[which(Choma4c$clust == clust)] <- xx$time_atpoint[jj]
  }
  Choma4b$clust[which(Choma4b$partid == part)] <- as.numeric(as.character(Choma4c$clust))
  Choma4b$points_inclust[which(Choma4b$partid == part)] <- Choma4c$points_inclust
  Choma4b$time_inclust[which(Choma4b$partid == part)] <- Choma4c$time_inclust
}
Choma4b$time_inclust_hour <- Choma4b$time_inclust*24
Choma4b$time_inclust_min <- Choma4b$time_inclust*24*60

### remove points greater than 1000 meters from center of cluster ####
Choma4b_short <- Choma4b[which(Choma4b$clust != 0),]
parts <- levels(as.factor(Choma4b_short$partid))
for(ii in 1:nlevels(as.factor(Choma4b_short$partid))){
  part <- parts[ii]
  temp <- Choma4b_short[which(Choma4b_short$partid == part),]
  out <- tapply(1:nrow(temp), as.factor(temp$clust), function(x) max(distm(as.matrix(temp[x,10:9], fun=distHaversine)))) ## meters
  clust_cut <- which(out > 2000)
  if(length(clust_cut) > 0){
    cluster_centers <- aggregate(temp[,10:9], list(Clust = as.factor(temp$clust)), median)
    for(jj in 1:length(clust_cut)){
      clust <- clust_cut[jj]
      temp2 <- temp[which(temp$clust == clust),]
      out2 <- distm(as.matrix(temp2[,10:9]), cluster_centers[which(cluster_centers$Clust==clust),2:3], fun=distHaversine)
      cuts <- which(out2 > 1000)
      Choma4b_short$clust[c(c(which(Choma4b_short$partid == part & Choma4b_short$clust == clust))[cuts])] <- 0
      Choma4b$clust[c(c(which(Choma4b$partid == part & Choma4b$clust == clust))[cuts])] <- 0
     
      Choma4b_short$points_inclust[c(which(Choma4b_short$partid == part & Choma4b_short$clust == clust))] <- nrow(Choma4b_short[c(which(Choma4b_short$partid == part & Choma4b_short$clust == clust)),])
      Choma4b$points_inclust[c(which(Choma4b$partid == part & Choma4b$clust == clust))] <- nrow(Choma4b[c(which(Choma4b$partid == part & Choma4b$clust == clust)),])
    }
  }
}
Choma4b_short2 <- Choma4b_short[which(Choma4b_short$clust != 0),]
Choma_short <- Choma4b_short2[!duplicated(Choma4b_short2[,c(1,21)]),]
##
###### limit number of clusters based on time spent and points ####
## by time spent
Choma_short_time <- Choma_short[order(Choma_short$partid, -Choma_short$time_inclust),]
### cut at 20 locations per person
Choma_short_time$in_out <- factor(NA, levels=c(0,1))
for(ii in 1:length(parts)){
  part <- parts[ii]
  ind <- which(Choma_short_time$partid == part)
  Choma_short_time$in_out[ind[1:20]] <- 1
  if(length(ind) > 20){
    Choma_short_time$in_out[ind[21:length(ind)]] <- 0
  }
}

## by number of points
Choma_short_points <- Choma_short[order(Choma_short$partid, -Choma_short$points_inclust),]
Choma_short_points$in_out <- factor(NA, levels=c(0,1))
for(ii in 1:length(parts)){
  part <- parts[ii]
  ind <- which(Choma_short_points$partid == part)
  Choma_short_points$in_out[ind[1:20]] <- 1
  if(length(ind) > 20){
    Choma_short_points$in_out[ind[21:length(ind)]] <- 0
  }
}
##
parts <- levels(as.factor(Choma_short$partid))
counts <- numeric(length=length(parts))
for(ii in 1:length(parts)){
  part <- parts[ii]
  counts[ii] <- length(sort(unique(c(Choma_short_points$clust[which(Choma_short_points$partid == part & Choma_short_points$in_out == 1)],
                                     Choma_short_time$clust[which(Choma_short_time$partid == part & Choma_short_time$in_out == 1)]))))
}
summary(as.factor(counts))
#  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
# 22 21 13 19 19 11 10 12  8 11  6  3  8 11 19 10  2  1  2 

Choma_short_points2 <- Choma_short_points[which(Choma_short_points$in_out == 1),]
Choma_short_time2 <- Choma_short_time[which(Choma_short_time$in_out == 1),]

##### remove clusters with less than 30 minutes or less than 5 points #####
Choma_short2 <- Choma_short
Choma_short2$keep <- factor(0, levels=c(0,1))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- sort(unique(c(Choma_short_points$clust[which(Choma_short_points$partid == part & Choma_short_points$in_out == 1)],
                          Choma_short_time$clust[which(Choma_short_time$partid == part & Choma_short_time$in_out == 1)])))
  Choma_short2$keep[which(Choma_short2$partid == part & Choma_short2$clust %in% clusts)] <- 1
}

Choma_short2b <- Choma_short2[which(Choma_short2$keep == 1),]

length(unique(Choma_short2b$partid[which(Choma_short2b$time_inclust_min < 20 | Choma_short2b$points_inclust < 5)])) ## 0
Choma_short2c <- Choma_short2b
a <- length(unique(Choma_short2c$partid[which(Choma_short2c$time_inclust_min < 30 | Choma_short2c$points_inclust < 5)])) ## 1
if(a == 0){
  Choma_short2d <- Choma_short2c}
if(a > 0){
  Choma_short2d <- Choma_short2c[-c(which(Choma_short2c$time_inclust_min < 30 | Choma_short2c$points_inclust < 5)),]
}
favstats(favstats(Choma_short2d$clust ~ Choma_short2d$partid)$n)
# min Q1 median Q3 max     mean      sd   n missing
#    1 8.25     14 21  25 14.16129 6.529111 62       0

#####
Choma5 <- Choma4b
Choma5$keep <- factor(0, levels=c(0,1))
parts <- levels(as.factor(Choma5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- Choma_short2d$clust[which(Choma_short2d$partid == part)]
  Choma5$keep[which(Choma5$partid == part & Choma5$clust %in% clusts)] <- 1
}

### define home cluster by distance from GPS of home listed #####
Choma_short2d$home_clust <- factor(0, levels=c(0,1))
Choma5$home_clust <- factor(0, levels=c(0,1))

parts <- levels(as.factor(Choma5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  min_dist <- min(Choma_short2d$hhdist_m_haversine[which(Choma_short2d$partid == part)])
  home_clust <- Choma_short2d$clust[which(Choma_short2d$partid == part & Choma_short2d$hhdist_m_haversine == min_dist)]
  Choma_short2d$home_clust[which(Choma_short2d$partid == part & Choma_short2d$clust == home_clust)] <- 1
  Choma5$home_clust[which(Choma5$partid == part & Choma5$clust == home_clust)] <- 1
}
#
##### define trips to locations #####

Choma5$time_atpoint_hour <- Choma5$time_atpoint*24
Choma5$time_atpoint_min <- Choma5$time_atpoint*24*60

Choma5$trip_count_new <- numeric(length=nrow(Choma5))
Choma5$trip_time_new <- numeric(length=nrow(Choma5))
parts <- levels(as.factor(Choma5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Choma5$clust[which(Choma5$partid == part & Choma5$keep == 1)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    inds <- which(Choma5$partid == part & Choma5$clust == clust)
    counts <- cumsum(c(1,diff(inds)>6))
    
    counts[which(counts %in% c(which(summary(as.factor(counts), maxsum=500) < 7)))]<- NA
    counts <- as.factor(counts)
    levels(counts) <- c(1:nlevels(counts))
    Choma5$trip_count_new[which(Choma5$partid == part & Choma5$clust == clust)] <- as.numeric(as.character(counts))
    Choma5c <- Choma5[which(Choma5$partid == part & Choma5$clust == clust),]
    if(length(which(!(is.na(Choma5c$trip_count_new))))){
      xx <- aggregate(time_atpoint ~ trip_count_new, Choma5c, sum)
      for(kk in 1:nrow(xx)){
        Choma5$trip_time_new[which(Choma5$partid == part & Choma5$clust == clust & Choma5$trip_count_new == kk)] <- xx$time_atpoint[kk]
      }
    }
  }
}
Choma5$trip_time_new[which(is.na(Choma5$trip_count_new))] <- NA
Choma5$trip_time_new_hour <- Choma5$trip_time_new*24
Choma5$trip_time_new_min <- Choma5$trip_time_new*24*60
#beep()
##### redfeine home cluster ###### 
### do 5 closest distances --> only if less than 1-2 km --> which spend most time at --> make that home cluster
Choma5$home_clust2 <- factor(0, levels=c(0,1))

parts <- levels(as.factor(Choma5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Choma5aa <- Choma5[which(Choma5$partid==part),]
  Choma5aa <- Choma5aa[order(Choma5aa$hhdist_m_haversine),]
  Choma5aa2 <- Choma5aa[!duplicated(Choma5aa$clust),]
  min_dists <- sort(Choma5aa2$hhdist_m_haversine[which(Choma5aa2$clust != 0 & Choma5aa2$keep==1)])
  min_dists2 <- min_dists[which(min_dists < 500)]
  clust_temp <- Choma5[which(Choma5$partid == part & Choma5$hhdist_m_haversine %in% min_dists2),]
  home_clust2 <- clust_temp$clust[which.max(clust_temp$time_inclust)]
  Choma5$home_clust2[which(Choma5$partid == part & Choma5$clust == home_clust2)] <- 1
}

### add one point per hour in location if <9 hrs missing overnight or <4hrs missing during the same day and points before and after gap in same cluster ####
Choma5$date_time_new <- as.POSIXct((Choma5$log_date_time_num*24*60*60), tz="Africa/Lusaka", origin="1899-12-30")
dates <- c(unlist(strsplit(as.character(Choma5$date_time_new)," "))[c(seq(from=1,to=c((nrow(Choma5)*2)-1), by=2))])
times <- c(unlist(strsplit(as.character(Choma5$date_time_new)," "))[c(seq(from=2,to=c(nrow(Choma5)*2), by=2))])
hours <- c(unlist(strsplit(times,":"))[c(seq(from=1,to=c((nrow(Choma5)*3)-2), by=3))])
Choma5$date_new <- dates
Choma5$time_new <- times
Choma5$time_new_hour <- as.numeric(as.character(hours))

Choma5_temp <- Choma5[order(Choma5$partid, Choma5$date_time_new),]

times <- 0
times2 <- 0
times_all <- 0

parts <- levels(as.factor(Choma5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Choma5_temp$clust[which(Choma5_temp$partid == part & Choma5_temp$keep == 1)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    temp <- Choma5_temp[which(Choma5_temp$partid == part & Choma5_temp$clust == clust & Choma5_temp$keep == 1),c(36,38)]
    for(kk in 1:(nrow(temp)-1)){
      if((as.Date(temp$date_new[kk])+1) == temp$date_new[kk + 1]){
        if(temp$time_new_hour[kk] != 23 | temp$time_new_hour[kk + 1] != 0){
          if((as.numeric(rownames(temp))[kk] + 1) == as.numeric(rownames(temp))[kk + 1]){
            length(times) <- 0
            length(times2) <- 0
            length(times_all) <- 0
            if((temp$time_new_hour[kk] != 23) & (temp$time_new_hour[kk + 1] != 0)){
              times <- seq((temp$time_new_hour[kk]+1),23,by=1)
              times2 <- seq(0,(temp$time_new_hour[kk + 1] - 1),by=1)
              times_all <- c(times,times2)
              add <- Choma5_temp[which(Choma5_temp$partid == part & Choma5_temp$clust == clust  & Choma5_temp$keep == 1)[kk],]
              add[,c(4:8,11:13,22:25,28:29,31:33,35)] <- NA
              adds <- as.data.frame(lapply(add, rep, length(times_all)))
              adds$date_new <- c(rep(temp$date_new[kk], times=length(times)), rep(temp$date_new[kk + 1], times=length(times2)))
              adds$time_new_hour <- times_all
              adds$time_new <- paste0(times_all, ":00:00")
            }
            if((temp$time_new_hour[kk] != 23) & (temp$time_new_hour[kk + 1] == 0)){
              times <- seq((temp$time_new_hour[kk]+1),23,by=1)
              times_all <- times
              add <- Choma5_temp[which(Choma5_temp$partid == part & Choma5_temp$clust == clust  & Choma5_temp$keep == 1)[kk],]
              add[,c(4:8,11:13,22:25,28:29,31:33,35)] <- NA
              adds <- as.data.frame(lapply(add, rep, length(times_all)))
              adds$date_new <- c(rep(temp$date_new[kk], times=length(times)))
              adds$time_new_hour <- times_all
              adds$time_new <- paste0(times_all, ":00:00")
            }
            if((temp$time_new_hour[kk] == 23) & (temp$time_new_hour[kk + 1] != 0)){
              times2 <- seq(0,(temp$time_new_hour[kk + 1] - 1),by=1)
              times_all <- times2
              add <- Choma5_temp[which(Choma5_temp$partid == part & Choma5_temp$clust == clust  & Choma5_temp$keep == 1)[kk],]
              add[,c(4:8,11:13,22:25,28:29,31:33,35)] <- NA
              adds <- as.data.frame(lapply(add, rep, length(times_all)))
              adds$date_new <- c(rep(temp$date_new[kk + 1], times=length(times2)))
              adds$time_new_hour <- times_all
              adds$time_new <- paste0(times_all, ":00:00")
            }
            if(length(times_all) < 9){
              Choma5_temp <- rbind(Choma5_temp,adds)
            }
          }
        }
      }
      if((temp$date_new[kk] == temp$date_new[kk + 1]) & (temp$time_new_hour[kk] != temp$time_new_hour[kk + 1])){
        if((temp$time_new_hour[kk] + 1) != temp$time_new_hour[kk + 1]){
          if((as.numeric(rownames(temp))[kk] + 1) == as.numeric(rownames(temp))[kk + 1]){
            add <- Choma5_temp[which(Choma5_temp$partid == part & Choma5_temp$clust == clust  & Choma5_temp$keep == 1)[kk],]
            add[,c(4:8,11:13,22:25,28:29,31:33,35)] <- NA
            count <- (temp$time_new_hour[kk + 1] - temp$time_new_hour[kk])-1
            if(count < 4){
              count2 <- seq((temp$time_new_hour[kk] + 1), by=1, length.out=count)
              adds <- as.data.frame(lapply(add, rep, count))
              adds$time_new_hour <- count2
              adds$time_new <- paste0(count2, ":00:00")
              Choma5_temp <- rbind(Choma5_temp,adds)
            }
          }
        }
      }
    }
  }
}
Choma5_temp$time_of_day <- character(length=nrow(Choma5_temp))
Choma5_temp$time_of_day[which(Choma5_temp$time_new_hour >= 6 & Choma5_temp$time_new_hour <= 11)] <- "morning"  ## 6am-11am
Choma5_temp$time_of_day[which(Choma5_temp$time_new_hour >= 12 & Choma5_temp$time_new_hour <= 17)] <- "midday"  ## 12pm-5pm
Choma5_temp$time_of_day[which(Choma5_temp$time_new_hour >= 18 & Choma5_temp$time_new_hour <= 23)] <- "evening"  ## 6pm-11pm
Choma5_temp$time_of_day[which(Choma5_temp$time_new_hour <= 5)] <- "night"  ## 12am-5am
Choma5_temp$time_of_day <- factor(Choma5_temp$time_of_day, levels=c("morning","midday","evening","night"))

Choma5b <- droplevels(Choma5_temp[which(Choma5_temp$keep == 1 & !(is.na(Choma5_temp$trip_count_new))),])

##
###### calculate new values for time per trip and per cluster #####
Choma5b2 <- Choma5b[order(Choma5b$partid, Choma5b$date_new, Choma5b$time_new_hour, Choma5b$time_new),]

Choma5b2$time_atpoint_new <- as.numeric(NA)
Choma5b2$trip_time_new2 <- as.numeric(NA)
Choma5b2$trip_time_new2_morning <- as.numeric(NA)
Choma5b2$trip_time_new2_midday <- as.numeric(NA)
Choma5b2$trip_time_new2_evening <- as.numeric(NA)
Choma5b2$trip_time_new2_night <- as.numeric(NA)
Choma5b2$clust_trips_time_new <- as.numeric(NA)
Choma5b2$clust_trips_time_new_morning <- as.numeric(NA)
Choma5b2$clust_trips_time_new_midday <- as.numeric(NA)
Choma5b2$clust_trips_time_new_evening <- as.numeric(NA)
Choma5b2$clust_trips_time_new_night <- as.numeric(NA)

parts <- levels(as.factor(Choma5b2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Choma5b2$clust[which(Choma5b2$partid == part)]))
  for(jj in 1:length(clusts)){
    clust <- clusts[jj]
    trips <- levels(as.factor(Choma5b2$trip_count_new[which(Choma5b2$partid == part & Choma5b2$clust == clust)]))
    for(kk in 1:length(trips)){
      trip <- trips[kk]
      temp2 <- Choma5b2[which(Choma5b2$partid == part & Choma5b2$clust == clust & Choma5b2$trip_count_new == trip),]
      time_at_points_middle <- diff(c(as.POSIXct(paste(temp2$date_new, temp2$time_new), format="%Y-%m-%d %H:%M:%S", tz="Africa/Lusaka")),2)
      if(units(time_at_points_middle) == "hours"){
        time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/24
      }
      if(units(time_at_points_middle) == "mins"){
        time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/60/24
      }
      if(units(time_at_points_middle) == "secs"){
        time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/60/60/24
      }
      time_at_points <- c(as.numeric(difftime(as.POSIXct(paste(temp2$date_new[2], temp2$time_new[2])), as.POSIXct(paste(temp2$date_new[1], temp2$time_new[1])), unit="min"))/2/60/24, 
                          time_at_points_middle2,
                          as.numeric(difftime(as.POSIXct(paste(temp2$date_new[c(nrow(temp2))], temp2$time_new[c(nrow(temp2))])), 
                                              as.POSIXct(paste(temp2$date_new[c(nrow(temp2)-1)], temp2$time_new[c(nrow(temp2)-1)])), 
                                              unit="min"))/2/60/24)
      Choma5b2$time_atpoint_new[which(Choma5b2$partid == part & Choma5b2$clust == clust & Choma5b2$trip_count_new == trip)] <- time_at_points
    }
    temp3 <- Choma5b2[which(Choma5b2$partid == part & Choma5b2$clust == clust),]
    xx <- aggregate(time_atpoint_new ~ trip_count_new, temp3, sum)
    for(kk in 1:nrow(xx)){
      Choma5b2$trip_time_new2[which(Choma5b2$partid == part & Choma5b2$clust == clust & Choma5b2$trip_count_new == kk)] <- xx$time_atpoint[kk]
    }
    xx2 <- aggregate(time_atpoint_new ~ trip_count_new + time_of_day, temp3, sum)
    xx2_morning <- xx2[which(xx2$time_of_day == "morning"),]
    xx2_midday <- xx2[which(xx2$time_of_day == "midday"),]
    xx2_evening <- xx2[which(xx2$time_of_day == "evening"),]
    xx2_night <- xx2[which(xx2$time_of_day == "night"),]
    for(kk in 1:nrow(xx2_morning)){
      Choma5b2$trip_time_new2_morning[which(Choma5b2$partid == part & Choma5b2$clust == clust & Choma5b2$trip_count_new == xx2_morning$trip_count_new[kk])] <- xx2_morning$time_atpoint_new[kk]
    }
    for(kk in 1:nrow(xx2_midday)){
      Choma5b2$trip_time_new2_midday[which(Choma5b2$partid == part & Choma5b2$clust == clust & Choma5b2$trip_count_new == xx2_midday$trip_count_new[kk])] <- xx2_midday$time_atpoint_new[kk]
    }
    for(kk in 1:nrow(xx2_evening)){
      Choma5b2$trip_time_new2_evening[which(Choma5b2$partid == part & Choma5b2$clust == clust & Choma5b2$trip_count_new == xx2_evening$trip_count_new[kk])] <- xx2_evening$time_atpoint_new[kk]
    }
    for(kk in 1:nrow(xx2_night)){
      Choma5b2$trip_time_new2_night[which(Choma5b2$partid == part & Choma5b2$clust == clust & Choma5b2$trip_count_new == xx2_night$trip_count_new[kk])] <- xx2_night$time_atpoint_new[kk]
    }
  }
  temp4 <- Choma5b2[which(Choma5b2$partid == part),]
  temp4b <- temp4[!duplicated(temp4[c(1,21,30)]),]
  xx <- aggregate(trip_time_new2 ~ clust, temp4b, sum)
  for(kk in 1:nrow(xx)){
    Choma5b2$clust_trips_time_new[which(Choma5b2$partid == part & Choma5b2$clust == xx$clust[kk] & !(is.na(Choma5b2$trip_count_new)))] <- xx$trip_time_new2[kk]
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_morning))),])>0){
    morning <- aggregate(trip_time_new2_morning ~ clust, temp4b, sum)
    if(nrow(morning) > 0){
      for(kk in 1:nrow(morning)){
        Choma5b2$clust_trips_time_new_morning[which(Choma5b2$partid == part & Choma5b2$clust == morning$clust[kk] & !(is.na(Choma5b2$trip_count_new)))] <- morning$trip_time_new2_morning[kk]
      }
    }
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_midday))),])>0){
    midday <- aggregate(trip_time_new2_midday ~ clust, temp4b, sum)
    if(nrow(midday) > 0){
      for(kk in 1:nrow(midday)){
        Choma5b2$clust_trips_time_new_midday[which(Choma5b2$partid == part & Choma5b2$clust == midday$clust[kk] & !(is.na(Choma5b2$trip_count_new)))] <- midday$trip_time_new2_midday[kk]
      }
    }
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_evening))),])>0){
    evening <- aggregate(trip_time_new2_evening ~ clust, temp4b, sum)
    if(nrow(evening) > 0){
      for(kk in 1:nrow(evening)){
        Choma5b2$clust_trips_time_new_evening[which(Choma5b2$partid == part & Choma5b2$clust == evening$clust[kk] & !(is.na(Choma5b2$trip_count_new)))] <- evening$trip_time_new2_evening[kk]
      }
    }
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_night))),])>0){
    night <- aggregate(trip_time_new2_night ~ clust, temp4b, sum)
    if(nrow(night) > 0){
      for(kk in 1:nrow(night)){
        Choma5b2$clust_trips_time_new_night[which(Choma5b2$partid == part & Choma5b2$clust == night$clust[kk] & !(is.na(Choma5b2$trip_count_new)))] <- night$trip_time_new2_night[kk]
      }
    }
  }
}


Choma5b2a <- Choma5b2[order(Choma5b2$partid, Choma5b2$clust, Choma5b2$trip_count_new),]
Choma5b2a2 <- Choma5b2a[which(!(is.na(Choma5b2a$trip_time_new2))),]


Choma5b2a2$part_trips_time_new <- as.numeric(0)
Choma5b2a2$part_trips_time_new_morning <- as.numeric(0)
Choma5b2a2$part_trips_time_new_midday <- as.numeric(0)
Choma5b2a2$part_trips_time_new_evening <- as.numeric(0)
Choma5b2a2$part_trips_time_new_night <- as.numeric(0)
temp <- Choma5b2a2[!duplicated(Choma5b2a2[c(1,21)]),]
tot <- aggregate(clust_trips_time_new ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Choma5b2a2$part_trips_time_new[which(Choma5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new[ii]
}
tot <- aggregate(clust_trips_time_new_morning ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Choma5b2a2$part_trips_time_new_morning[which(Choma5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_morning[ii]
}
tot <- aggregate(clust_trips_time_new_midday ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Choma5b2a2$part_trips_time_new_midday[which(Choma5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_midday[ii]
}
tot <- aggregate(clust_trips_time_new_evening ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Choma5b2a2$part_trips_time_new_evening[which(Choma5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_evening[ii]
}
tot <- aggregate(clust_trips_time_new_night ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Choma5b2a2$part_trips_time_new_night[which(Choma5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_night[ii]
}

#######
Choma5_new <- Choma5b2a2[,c(1:3,14,15,4,12,5:11,29,16:20,27,21:23,25,30,31,33,34,35:55,26)]

### take out those with <15 minutes in cluster or <15 minute trips that may be here now because of clusters that were >1000 meter radius
Choma5_new2 <- Choma5_new[which(Choma5_new$clust_trips_time_new > 0.01),]
Choma5_new2b <- Choma5_new2[which(Choma5_new2$trip_time_new2 > 0.01),]

Choma5_clusts_shorter <- Choma5_new2b[order(Choma5_new2b$partid, Choma5_new2b$clust),]
Choma5_clusts_shorter2 <- Choma5_clusts_shorter[!duplicated(Choma5_clusts_shorter[,c(1,22)]),c(1:7,16,17,12,13,19:25,29,30:35,46:50,41:45,51)]

Choma5_clust_trips_shorter <- Choma5_new2b[order(Choma5_new2b$partid, Choma5_new2b$clust, Choma5_new2b$trip_count_new),]
Choma5_clust_trips_shorter2 <- Choma5_clust_trips_shorter[!duplicated(Choma5_clust_trips_shorter[,c(1,22,26)]),c(1:7,16,17,12,13,19:22,24,26:28,29,30:35,46:50,41:45,36:40,51)]
##
#####
setwd()
write.csv(Choma4, "Choma4.csv", row.names = FALSE)
write.csv(Choma4b, "Choma4b.csv", row.names = FALSE)
write.csv(Choma5, "Choma5.csv", row.names = FALSE)
write.csv(Choma5_new2b, "Choma5_new.csv", row.names = FALSE)

write.csv(Choma5_clusts_shorter2, "Choma5_clusts_shorter2.csv", row.names = FALSE)
write.csv(Choma5_clust_trips_shorter2, "Choma5_clust_trips_shorter2.csv", row.names = FALSE)
##

############ 
############ 
############ Nchelenge ###############
############ 
############ 
##### read in data ####
setwd()
Nchelenge4 <- read.csv("Nchelenge4.csv", header=TRUE, sep=",")
##### remove those with less than 500 points #####
Nchelenge4b <- Nchelenge4[,-c(19:20)]
Nchelenge4b <- droplevels(Nchelenge4b[which(Nchelenge4b$tot_points > 500),])
Nchelenge4b <- Nchelenge4b[order(Nchelenge4b$partid, Nchelenge4b$point_num),]
#
### fix those with excessive points taken affecting clusters (points every 10 seconds or less when 15+ in a row)--> make once a min #####
parts_too_close <- as.factor(names(which(table(as.factor(Nchelenge4b$partid[which((Nchelenge4b$time_atpoint) < 0.00011574)])) > 10)))
for(ii in 1:length(parts_too_close)){
  ind_keep <- NA
  ind_rem <- NA
  part <- parts_too_close[ii]
  temp <- Nchelenge4b[which(Nchelenge4b$partid == part),]
  inds <- which(diff(temp$log_date_time_num) < 0.00011574)
  counts <- cumsum(c(1,diff(inds)!=1))
  counts[which(counts %in% c(which(summary(as.factor(counts), maxsum=1000) < 15)))]<- NA
  counts <- as.factor(counts)
  for(jj in 1:nlevels(counts)){
    count <- levels(counts)[jj]
    inds2 <- inds[which(counts == count)]
    ind_keep <- c(ind_keep,inds2[seq(from=1,to=length(inds2), by=12)])
    ind_rem <- c(ind_rem,inds2[-c(seq(from=1,to=length(inds2), by=12))])
  }
  ind_rem <- ind_rem[-1]
  Nchelenge4b <- Nchelenge4b[-c(which(Nchelenge4b$partid == part)[ind_rem]),]
  temp <- Nchelenge4b[which(Nchelenge4b$partid == part),]
  time_at_points <- c(c(diff(c(temp$log_date_time_num[1], temp$log_date_time_num[2]))/2),
                      c(diff(c(temp$log_date_time_num),2)/2),
                      c(diff(c(temp$log_date_time_num[c(nrow(temp)-1)], temp$log_date_time_num[c(nrow(temp))]))/2))
  Nchelenge4b$time_atpoint[which(Nchelenge4b$partid == part)] <- time_at_points
  Nchelenge4b$total_time[which(Nchelenge4b$partid == part)] <- sum(time_at_points)
}
#
#### define clusters with HDBSCAN and DBSCAN #####
Nchelenge4b$clust <- numeric(length=nrow(Nchelenge4b))
Nchelenge4b$points_inclust <- numeric(length=nrow(Nchelenge4b))
Nchelenge4b$time_inclust <- numeric(length=nrow(Nchelenge4b))

parts <- levels(as.factor(Nchelenge4b$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Nchelenge4c <- Nchelenge4b[which(Nchelenge4b$partid == part),]
  xy <- SpatialPointsDataFrame(
    matrix(c(Nchelenge4c$log_long, Nchelenge4c$log_lat), ncol=2), data.frame(ID=seq(1:nrow(Nchelenge4c))),
    proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  mdist <- distm(xy)
  bb <- hdbscan(as.dist(mdist), minPts = 15) 
  ##  if more than 2 clusters (and noise), all good
  if(length(unique(bb$cluster)) > 3){
    xy$clust <- bb$cluster
  }
  ## if 2 clusters (and noise)
  else if(length(unique(bb$cluster)) == 3){
    res <- .Call('_dbscan_computeStability', bb$hc, 15, compute_glosh=TRUE)
    if(res$'0'$stability < (res$'1'$stability + res$'2'$stability)){ ## 2 clusters are good, keep going as normal
      xy$clust <- bb$cluster
    }
    else{ ## root is better --> still want to find outliers though
      res2 <- dbscan(as.dist(mdist), eps = 25, minPts = 6)
      clust_temp <- as.numeric(as.character(as.factor(res2$cluster)[(which.max(summary(as.factor(res2$cluster))))])) ## not actual cluster # with which.max since minPts > 1
      xy$clust <- 0
      xy$clust[which(res2$cluster == clust_temp)] <- 1
    }
  }
  else if(length(unique(bb$cluster)) == 2){ ## either 1 clust and outliers or 2 clusts and no outliers; either way, should be good
    xy$clust <- bb$cluster
  }
  else if(length(unique(bb$cluster)) == 1){ ## only 1 clust found but all considered 'outliers'
    res2 <- dbscan(as.dist(mdist), eps = 25, minPts = 1)
    clust_temp <- which.max(summary(as.factor(res2$cluster)))
    xy$clust <- 0
    xy$clust[which(res2$cluster == clust_temp)] <- 1
  }
  if(length(which(bb$cluster == 0))/length(bb$cluster) > 0.5){
    res2 <- dbscan(as.dist(mdist), eps = 25, minPts = 1)
    clust_temp <- which.max(summary(as.factor(res2$cluster))) ## actual clust # since no noise points (since minPts =1)
    xy$clust <- 0
    xy$clust[which(res2$cluster == clust_temp)] <- 1
  }
  xy2 <- as.data.frame(xy)
  Nchelenge4c$clust <- factor(as.numeric(xy2$clust))
  xx <- aggregate(time_atpoint ~ clust, Nchelenge4c, sum)
  for(jj in 1:nrow(xx)){
    clust <- xx$clust[jj]
    
    Nchelenge4c$points_inclust[which(Nchelenge4c$clust == clust)] <- length(which(xy2$clust == clust))
    Nchelenge4c$time_inclust[which(Nchelenge4c$clust == clust)] <- xx$time_atpoint[jj]
  }
  Nchelenge4b$clust[which(Nchelenge4b$partid == part)] <- as.numeric(as.character(Nchelenge4c$clust))
  Nchelenge4b$points_inclust[which(Nchelenge4b$partid == part)] <- Nchelenge4c$points_inclust
  Nchelenge4b$time_inclust[which(Nchelenge4b$partid == part)] <- Nchelenge4c$time_inclust
}
Nchelenge4b$time_inclust_hour <- Nchelenge4b$time_inclust*24
Nchelenge4b$time_inclust_min <- Nchelenge4b$time_inclust*24*60
Nchelenge_short1 <- Nchelenge4b[!duplicated(Nchelenge4b[,c(1,21)]),]

Nchelenge_short <- Nchelenge_short1[which(Nchelenge_short1$clust != 0),]
#beep()

### remove points greater than 1000 meters from center of cluster ####
Nchelenge4b_short <- Nchelenge4b[which(Nchelenge4b$clust != 0),]
parts <- levels(as.factor(Nchelenge4b_short$partid))
for(ii in 1:nlevels(as.factor(Nchelenge4b_short$partid))){
  part <- parts[ii]
  temp <- Nchelenge4b_short[which(Nchelenge4b_short$partid == part),]
  out <- tapply(1:nrow(temp), as.factor(temp$clust), function(x) max(distm(as.matrix(temp[x,10:9], fun=distHaversine)))) ## meters
  clust_cut <- which(out > 2000)
  if(length(clust_cut) > 0){
    cluster_centers <- aggregate(temp[,10:9], list(Clust = as.factor(temp$clust)), median)
    for(jj in 1:length(clust_cut)){
      clust <- clust_cut[jj]
      temp2 <- temp[which(temp$clust == clust),]
      out2 <- distm(as.matrix(temp2[,10:9]), cluster_centers[which(cluster_centers$Clust==clust),2:3], fun=distHaversine)
      cuts <- which(out2 > 1000)
      Nchelenge4b_short$clust[c(c(which(Nchelenge4b_short$partid == part & Nchelenge4b_short$clust == clust))[cuts])] <- 0
      Nchelenge4b$clust[c(c(which(Nchelenge4b$partid == part & Nchelenge4b$clust == clust))[cuts])] <- 0
      
      Nchelenge4b_short$points_inclust[c(which(Nchelenge4b_short$partid == part & Nchelenge4b_short$clust == clust))] <- nrow(Nchelenge4b_short[c(which(Nchelenge4b_short$partid == part & Nchelenge4b_short$clust == clust)),])
      Nchelenge4b$points_inclust[c(which(Nchelenge4b$partid == part & Nchelenge4b$clust == clust))] <- nrow(Nchelenge4b[c(which(Nchelenge4b$partid == part & Nchelenge4b$clust == clust)),])
    }
  }
}
Nchelenge4b_short2 <- Nchelenge4b_short[which(Nchelenge4b_short$clust != 0),]
Nchelenge_short <- Nchelenge4b_short2[!duplicated(Nchelenge4b_short2[,c(1,21)]),]
##
###### limit number of clusters based on time spent and points ####
## by time spent
Nchelenge_short_time <- Nchelenge_short[order(Nchelenge_short$partid, -Nchelenge_short$time_inclust),]
### cut at 20 locations per person
Nchelenge_short_time$in_out <- factor(NA, levels=c(0,1))
for(ii in 1:length(parts)){
  part <- parts[ii]
  ind <- which(Nchelenge_short_time$partid == part)
  Nchelenge_short_time$in_out[ind[1:20]] <- 1
  if(length(ind) > 20){
    Nchelenge_short_time$in_out[ind[21:length(ind)]] <- 0
  }
}

## by number of points
Nchelenge_short_points <- Nchelenge_short[order(Nchelenge_short$partid, -Nchelenge_short$points_inclust),]
Nchelenge_short_points$in_out <- factor(NA, levels=c(0,1))
for(ii in 1:length(parts)){
  part <- parts[ii]
  ind <- which(Nchelenge_short_points$partid == part)
  Nchelenge_short_points$in_out[ind[1:20]] <- 1
  if(length(ind) > 20){
    Nchelenge_short_points$in_out[ind[21:length(ind)]] <- 0
  }
}
##
parts <- levels(as.factor(Nchelenge_short$partid))
counts <- numeric(length=length(parts))
for(ii in 1:length(parts)){
  part <- parts[ii]
  counts[ii] <- length(sort(unique(c(Nchelenge_short_points$clust[which(Nchelenge_short_points$partid == part & Nchelenge_short_points$in_out == 1)],
                                     Nchelenge_short_time$clust[which(Nchelenge_short_time$partid == part & Nchelenge_short_time$in_out == 1)]))))
}
summary(as.factor(counts))
#  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
# 22 21 13 19 19 11 10 12  8 11  6  3  8 11 19 10  2  1  2 

Nchelenge_short_points2 <- Nchelenge_short_points[which(Nchelenge_short_points$in_out == 1),]
Nchelenge_short_time2 <- Nchelenge_short_time[which(Nchelenge_short_time$in_out == 1),]

##### remove clusters with less than 30 minutes or less than 5 points #####
Nchelenge_short2 <- Nchelenge_short
Nchelenge_short2$keep <- factor(0, levels=c(0,1))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- sort(unique(c(Nchelenge_short_points$clust[which(Nchelenge_short_points$partid == part & Nchelenge_short_points$in_out == 1)],
                          Nchelenge_short_time$clust[which(Nchelenge_short_time$partid == part & Nchelenge_short_time$in_out == 1)])))
  Nchelenge_short2$keep[which(Nchelenge_short2$partid == part & Nchelenge_short2$clust %in% clusts)] <- 1
}

Nchelenge_short2b <- Nchelenge_short2[which(Nchelenge_short2$keep == 1),]

length(unique(Nchelenge_short2b$partid[which(Nchelenge_short2b$time_inclust_min < 20 | Nchelenge_short2b$points_inclust < 5)])) ## 0
Nchelenge_short2c <- Nchelenge_short2b
a <- length(unique(Nchelenge_short2c$partid[which(Nchelenge_short2c$time_inclust_min < 30 | Nchelenge_short2c$points_inclust < 5)])) ## 1
if(a == 0){
  Nchelenge_short2d <- Nchelenge_short2c}
if(a > 0){
  Nchelenge_short2d <- Nchelenge_short2c[-c(which(Nchelenge_short2c$time_inclust_min < 30 | Nchelenge_short2c$points_inclust < 5)),]
}
favstats(favstats(Nchelenge_short2d$clust ~ Nchelenge_short2d$partid)$n)
# min Q1 median Q3 max     mean      sd   n missing
#   2  4      7 14  20 8.600962 5.139541 208       0

#####
Nchelenge5 <- Nchelenge4b
Nchelenge5$keep <- factor(0, levels=c(0,1))
parts <- levels(as.factor(Nchelenge5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- Nchelenge_short2d$clust[which(Nchelenge_short2d$partid == part)]
  Nchelenge5$keep[which(Nchelenge5$partid == part & Nchelenge5$clust %in% clusts)] <- 1
}


### define home cluster by distance from GPS of home listed #####
Nchelenge_short2d$home_clust <- factor(0, levels=c(0,1))
Nchelenge5$home_clust <- factor(0, levels=c(0,1))

parts <- levels(as.factor(Nchelenge5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  min_dist <- min(Nchelenge_short2d$hhdist_m_haversine[which(Nchelenge_short2d$partid == part)])
  home_clust <- Nchelenge_short2d$clust[which(Nchelenge_short2d$partid == part & Nchelenge_short2d$hhdist_m_haversine == min_dist)]
  Nchelenge_short2d$home_clust[which(Nchelenge_short2d$partid == part & Nchelenge_short2d$clust == home_clust)] <- 1
  Nchelenge5$home_clust[which(Nchelenge5$partid == part & Nchelenge5$clust == home_clust)] <- 1
}
#
##### define trips to locations #####
Nchelenge5$time_atpoint_hour <- Nchelenge5$time_atpoint*24
Nchelenge5$time_atpoint_min <- Nchelenge5$time_atpoint*24*60

Nchelenge5$trip_count_new <- numeric(length=nrow(Nchelenge5))
Nchelenge5$trip_time_new <- numeric(length=nrow(Nchelenge5))
parts <- levels(as.factor(Nchelenge5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Nchelenge5$clust[which(Nchelenge5$partid == part & Nchelenge5$keep == 1)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    inds <- which(Nchelenge5$partid == part & Nchelenge5$clust == clust)
    counts <- cumsum(c(1,diff(inds)>6))
    
    counts[which(counts %in% c(which(summary(as.factor(counts), maxsum=500) < 7)))]<- NA
    counts <- as.factor(counts)
    levels(counts) <- c(1:nlevels(counts))
    Nchelenge5$trip_count_new[which(Nchelenge5$partid == part & Nchelenge5$clust == clust)] <- as.numeric(as.character(counts))
    Nchelenge5c <- Nchelenge5[which(Nchelenge5$partid == part & Nchelenge5$clust == clust),]
    if(length(which(!(is.na(Nchelenge5c$trip_count_new))))){
      xx <- aggregate(time_atpoint ~ trip_count_new, Nchelenge5c, sum)
      for(kk in 1:nrow(xx)){
        Nchelenge5$trip_time_new[which(Nchelenge5$partid == part & Nchelenge5$clust == clust & Nchelenge5$trip_count_new == kk)] <- xx$time_atpoint[kk]
      }
    }
  }
}
Nchelenge5$trip_time_new[which(is.na(Nchelenge5$trip_count_new))] <- NA
Nchelenge5$trip_time_new_hour <- Nchelenge5$trip_time_new*24
Nchelenge5$trip_time_new_min <- Nchelenge5$trip_time_new*24*60
#beep()
##### redfeine home cluster ###### 
### do 5 closest distances --> only if less than 1-2 km --> which spend most time at --> make that home cluster
Nchelenge5$home_clust2 <- factor(0, levels=c(0,1))

parts <- levels(as.factor(Nchelenge5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Nchelenge5aa <- Nchelenge5[which(Nchelenge5$partid==part),]
  Nchelenge5aa <- Nchelenge5aa[order(Nchelenge5aa$hhdist_m_haversine),]
  Nchelenge5aa2 <- Nchelenge5aa[!duplicated(Nchelenge5aa$clust),]
  min_dists <- sort(Nchelenge5aa2$hhdist_m_haversine[which(Nchelenge5aa2$clust != 0 & Nchelenge5aa2$keep==1)])
  min_dists2 <- min_dists[which(min_dists < 500)]
  clust_temp <- Nchelenge5[which(Nchelenge5$partid == part & Nchelenge5$hhdist_m_haversine %in% min_dists2),]
  home_clust2 <- clust_temp$clust[which.max(clust_temp$time_inclust)]
  Nchelenge5$home_clust2[which(Nchelenge5$partid == part & Nchelenge5$clust == home_clust2)] <- 1
}
#
### add one point per hour in location if <9 hrs missing overnight or <4hrs missing during the same day and points before and after gap in same cluster ####
Nchelenge5$date_time_new <- as.POSIXct((Nchelenge5$log_date_time_num*24*60*60), tz="Africa/Lusaka", origin="1899-12-30")
dates <- c(unlist(strsplit(as.character(Nchelenge5$date_time_new)," "))[c(seq(from=1,to=c((nrow(Nchelenge5)*2)-1), by=2))])
times <- c(unlist(strsplit(as.character(Nchelenge5$date_time_new)," "))[c(seq(from=2,to=c(nrow(Nchelenge5)*2), by=2))])
hours <- c(unlist(strsplit(times,":"))[c(seq(from=1,to=c((nrow(Nchelenge5)*3)-2), by=3))])
Nchelenge5$date_new <- dates
Nchelenge5$time_new <- times
Nchelenge5$time_new_hour <- as.numeric(as.character(hours))

Nchelenge5_temp <- Nchelenge5[order(Nchelenge5$partid, Nchelenge5$date_time_new),]
Nchelenge5_temp <- Nchelenge5_temp[,c(1:19,21:28,20,29:37)]

times <- 0
times2 <- 0
times_all <- 0

parts <- levels(as.factor(Nchelenge5$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Nchelenge5_temp$clust[which(Nchelenge5_temp$partid == part & Nchelenge5_temp$keep == 1)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    temp <- Nchelenge5_temp[which(Nchelenge5_temp$partid == part & Nchelenge5_temp$clust == clust & Nchelenge5_temp$keep == 1),c(35,37)]
    for(kk in 1:(nrow(temp)-1)){
      if((as.Date(temp$date_new[kk])+1) == temp$date_new[kk + 1]){
        if(temp$time_new_hour[kk] != 23 | temp$time_new_hour[kk + 1] != 0){
          if((as.numeric(rownames(temp))[kk] + 1) == as.numeric(rownames(temp))[kk + 1]){
            length(times) <- 0
            length(times2) <- 0
            length(times_all) <- 0
            if((temp$time_new_hour[kk] != 23) & (temp$time_new_hour[kk + 1] != 0)){
              times <- seq((temp$time_new_hour[kk]+1),23,by=1)
              times2 <- seq(0,(temp$time_new_hour[kk + 1] - 1),by=1)
              times_all <- c(times,times2)
              add <- Nchelenge5_temp[which(Nchelenge5_temp$partid == part & Nchelenge5_temp$clust == clust  & Nchelenge5_temp$keep == 1)[kk],]
              add[,c(4:8,11:13,21:24,27:28,30:32,34)] <- NA
              adds <- as.data.frame(lapply(add, rep, length(times_all)))
              adds$date_new <- c(rep(temp$date_new[kk], times=length(times)), rep(temp$date_new[kk + 1], times=length(times2)))
              adds$time_new_hour <- times_all
              adds$time_new <- paste0(times_all, ":00:00")
            }
            if((temp$time_new_hour[kk] != 23) & (temp$time_new_hour[kk + 1] == 0)){
              times <- seq((temp$time_new_hour[kk]+1),23,by=1)
              times_all <- times
              add <- Nchelenge5_temp[which(Nchelenge5_temp$partid == part & Nchelenge5_temp$clust == clust  & Nchelenge5_temp$keep == 1)[kk],]
              add[,c(4:8,11:13,21:24,27:28,30:32,34)] <- NA
              adds <- as.data.frame(lapply(add, rep, length(times_all)))
              adds$date_new <- c(rep(temp$date_new[kk], times=length(times)))
              adds$time_new_hour <- times_all
              adds$time_new <- paste0(times_all, ":00:00")
            }
            if((temp$time_new_hour[kk] == 23) & (temp$time_new_hour[kk + 1] != 0)){
              times2 <- seq(0,(temp$time_new_hour[kk + 1] - 1),by=1)
              times_all <- times2
              add <- Nchelenge5_temp[which(Nchelenge5_temp$partid == part & Nchelenge5_temp$clust == clust  & Nchelenge5_temp$keep == 1)[kk],]
              add[,c(4:8,11:13,21:24,27:28,30:32,34)] <- NA
              adds <- as.data.frame(lapply(add, rep, length(times_all)))
              adds$date_new <- c(rep(temp$date_new[kk + 1], times=length(times2)))
              adds$time_new_hour <- times_all
              adds$time_new <- paste0(times_all, ":00:00")
            }
            if(length(times_all) < 9){
              Nchelenge5_temp <- rbind(Nchelenge5_temp,adds)
            }
          }
        }
      }
      if((temp$date_new[kk] == temp$date_new[kk + 1]) & (temp$time_new_hour[kk] != temp$time_new_hour[kk + 1])){
        if((temp$time_new_hour[kk] + 1) != temp$time_new_hour[kk + 1]){
          if((as.numeric(rownames(temp))[kk] + 1) == as.numeric(rownames(temp))[kk + 1]){
            add <- Nchelenge5_temp[which(Nchelenge5_temp$partid == part & Nchelenge5_temp$clust == clust  & Nchelenge5_temp$keep == 1)[kk],]
            add[,c(4:8,11:13,21:24,27:28,30:32,34)] <- NA
            count <- (temp$time_new_hour[kk + 1] - temp$time_new_hour[kk])-1
            if(count < 4){
              count2 <- seq((temp$time_new_hour[kk] + 1), by=1, length.out=count)
              adds <- as.data.frame(lapply(add, rep, count))
              adds$time_new_hour <- count2
              adds$time_new <- paste0(count2, ":00:00")
              Nchelenge5_temp <- rbind(Nchelenge5_temp,adds)
            }
          }
        }
      }
    }
  }
}

Nchelenge5_temp$time_of_day <- character(length=nrow(Nchelenge5_temp))
Nchelenge5_temp$time_of_day[which(Nchelenge5_temp$time_new_hour >= 6 & Nchelenge5_temp$time_new_hour <= 11)] <- "morning"  ## 6am-11am
Nchelenge5_temp$time_of_day[which(Nchelenge5_temp$time_new_hour >= 12 & Nchelenge5_temp$time_new_hour <= 17)] <- "midday"  ## 12pm-5pm
Nchelenge5_temp$time_of_day[which(Nchelenge5_temp$time_new_hour >= 18 & Nchelenge5_temp$time_new_hour <= 23)] <- "evening"  ## 6pm-11pm
Nchelenge5_temp$time_of_day[which(Nchelenge5_temp$time_new_hour <= 5)] <- "night"  ## 12am-5am
Nchelenge5_temp$time_of_day <- factor(Nchelenge5_temp$time_of_day, levels=c("morning","midday","evening","night"))

Nchelenge5b <- droplevels(Nchelenge5_temp[which(Nchelenge5_temp$keep == 1 & !(is.na(Nchelenge5_temp$trip_count_new))),])

##
###### calculate new values for time per trip and per cluster #####

Nchelenge5b2 <- Nchelenge5b[order(Nchelenge5b$partid, Nchelenge5b$date_new, Nchelenge5b$time_new_hour, Nchelenge5b$time_new),]

Nchelenge5b2$time_atpoint_new <- as.numeric(NA)
Nchelenge5b2$trip_time_new2 <- as.numeric(NA)
Nchelenge5b2$trip_time_new2_morning <- as.numeric(NA)
Nchelenge5b2$trip_time_new2_midday <- as.numeric(NA)
Nchelenge5b2$trip_time_new2_evening <- as.numeric(NA)
Nchelenge5b2$trip_time_new2_night <- as.numeric(NA)
Nchelenge5b2$clust_trips_time_new <- as.numeric(NA)
Nchelenge5b2$clust_trips_time_new_morning <- as.numeric(NA)
Nchelenge5b2$clust_trips_time_new_midday <- as.numeric(NA)
Nchelenge5b2$clust_trips_time_new_evening <- as.numeric(NA)
Nchelenge5b2$clust_trips_time_new_night <- as.numeric(NA)

parts <- levels(as.factor(Nchelenge5b2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Nchelenge5b2$clust[which(Nchelenge5b2$partid == part)]))
  for(jj in 1:length(clusts)){
    clust <- clusts[jj]
    trips <- levels(as.factor(Nchelenge5b2$trip_count_new[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == clust)]))
    for(kk in 1:length(trips)){
      trip <- trips[kk]
      temp2 <- Nchelenge5b2[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == clust & Nchelenge5b2$trip_count_new == trip),]
      time_at_points_middle <- diff(c(as.POSIXct(paste(temp2$date_new, temp2$time_new), format="%Y-%m-%d %H:%M:%S", tz="Africa/Lusaka")),2)
      if(units(time_at_points_middle) == "hours"){
        time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/24
      }
      if(units(time_at_points_middle) == "mins"){
        time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/60/24
      }
      if(units(time_at_points_middle) == "secs"){
        time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/60/60/24
      }
      time_at_points <- c(as.numeric(difftime(as.POSIXct(paste(temp2$date_new[2], temp2$time_new[2])), as.POSIXct(paste(temp2$date_new[1], temp2$time_new[1])), unit="min"))/2/60/24, 
                          time_at_points_middle2,
                          as.numeric(difftime(as.POSIXct(paste(temp2$date_new[c(nrow(temp2))], temp2$time_new[c(nrow(temp2))])), 
                                              as.POSIXct(paste(temp2$date_new[c(nrow(temp2)-1)], temp2$time_new[c(nrow(temp2)-1)])), 
                                              unit="min"))/2/60/24)
      Nchelenge5b2$time_atpoint_new[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == clust & Nchelenge5b2$trip_count_new == trip)] <- time_at_points
    }
    temp3 <- Nchelenge5b2[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == clust),]
    xx <- aggregate(time_atpoint_new ~ trip_count_new, temp3, sum)
    for(kk in 1:nrow(xx)){
      Nchelenge5b2$trip_time_new2[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == clust & Nchelenge5b2$trip_count_new == kk)] <- xx$time_atpoint[kk]
    }
    xx2 <- aggregate(time_atpoint_new ~ trip_count_new + time_of_day, temp3, sum)
    xx2_morning <- xx2[which(xx2$time_of_day == "morning"),]
    xx2_midday <- xx2[which(xx2$time_of_day == "midday"),]
    xx2_evening <- xx2[which(xx2$time_of_day == "evening"),]
    xx2_night <- xx2[which(xx2$time_of_day == "night"),]
    for(kk in 1:nrow(xx2_morning)){
      Nchelenge5b2$trip_time_new2_morning[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == clust & Nchelenge5b2$trip_count_new == xx2_morning$trip_count_new[kk])] <- xx2_morning$time_atpoint_new[kk]
    }
    for(kk in 1:nrow(xx2_midday)){
      Nchelenge5b2$trip_time_new2_midday[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == clust & Nchelenge5b2$trip_count_new == xx2_midday$trip_count_new[kk])] <- xx2_midday$time_atpoint_new[kk]
    }
    for(kk in 1:nrow(xx2_evening)){
      Nchelenge5b2$trip_time_new2_evening[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == clust & Nchelenge5b2$trip_count_new == xx2_evening$trip_count_new[kk])] <- xx2_evening$time_atpoint_new[kk]
    }
    for(kk in 1:nrow(xx2_night)){
      Nchelenge5b2$trip_time_new2_night[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == clust & Nchelenge5b2$trip_count_new == xx2_night$trip_count_new[kk])] <- xx2_night$time_atpoint_new[kk]
    }
  }
  temp4 <- Nchelenge5b2[which(Nchelenge5b2$partid == part),]
  temp4b <- temp4[!duplicated(temp4[c(1,20,29)]),]
  xx <- aggregate(trip_time_new2 ~ clust, temp4b, sum)
  for(kk in 1:nrow(xx)){
    Nchelenge5b2$clust_trips_time_new[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == xx$clust[kk] & !(is.na(Nchelenge5b2$trip_count_new)))] <- xx$trip_time_new2[kk]
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_morning))),])>0){
    morning <- aggregate(trip_time_new2_morning ~ clust, temp4b, sum)
    if(nrow(morning) > 0){
      for(kk in 1:nrow(morning)){
        Nchelenge5b2$clust_trips_time_new_morning[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == morning$clust[kk] & !(is.na(Nchelenge5b2$trip_count_new)))] <- morning$trip_time_new2_morning[kk]
      }
    }
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_midday))),])>0){
    midday <- aggregate(trip_time_new2_midday ~ clust, temp4b, sum)
    if(nrow(midday) > 0){
      for(kk in 1:nrow(midday)){
        Nchelenge5b2$clust_trips_time_new_midday[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == midday$clust[kk] & !(is.na(Nchelenge5b2$trip_count_new)))] <- midday$trip_time_new2_midday[kk]
      }
    }
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_evening))),])>0){
    evening <- aggregate(trip_time_new2_evening ~ clust, temp4b, sum)
    if(nrow(evening) > 0){
      for(kk in 1:nrow(evening)){
        Nchelenge5b2$clust_trips_time_new_evening[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == evening$clust[kk] & !(is.na(Nchelenge5b2$trip_count_new)))] <- evening$trip_time_new2_evening[kk]
      }
    }
  }
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_night))),])>0){
    night <- aggregate(trip_time_new2_night ~ clust, temp4b, sum)
    if(nrow(night) > 0){
      for(kk in 1:nrow(night)){
        Nchelenge5b2$clust_trips_time_new_night[which(Nchelenge5b2$partid == part & Nchelenge5b2$clust == night$clust[kk] & !(is.na(Nchelenge5b2$trip_count_new)))] <- night$trip_time_new2_night[kk]
      }
    }
  }
}


Nchelenge5b2a <- Nchelenge5b2[order(Nchelenge5b2$partid, Nchelenge5b2$clust, Nchelenge5b2$trip_count_new),]
Nchelenge5b2a2 <- Nchelenge5b2a[which(!(is.na(Nchelenge5b2a$trip_time_new2))),]


Nchelenge5b2a2$part_trips_time_new <- as.numeric(0)
Nchelenge5b2a2$part_trips_time_new_morning <- as.numeric(0)
Nchelenge5b2a2$part_trips_time_new_midday <- as.numeric(0)
Nchelenge5b2a2$part_trips_time_new_evening <- as.numeric(0)
Nchelenge5b2a2$part_trips_time_new_night <- as.numeric(0)
temp <- Nchelenge5b2a2[!duplicated(Nchelenge5b2a2[c(1,20)]),]
tot <- aggregate(clust_trips_time_new ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Nchelenge5b2a2$part_trips_time_new[which(Nchelenge5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new[ii]
}
tot <- aggregate(clust_trips_time_new_morning ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Nchelenge5b2a2$part_trips_time_new_morning[which(Nchelenge5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_morning[ii]
}
tot <- aggregate(clust_trips_time_new_midday ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Nchelenge5b2a2$part_trips_time_new_midday[which(Nchelenge5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_midday[ii]
}
tot <- aggregate(clust_trips_time_new_evening ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Nchelenge5b2a2$part_trips_time_new_evening[which(Nchelenge5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_evening[ii]
}
tot <- aggregate(clust_trips_time_new_night ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Nchelenge5b2a2$part_trips_time_new_night[which(Nchelenge5b2a2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_night[ii]
}
##
##### AT PRIMARY OR SECONDARY HOME? #####
Nchelenge <- read.csv("Nchelenge_AllPoints.csv", header=TRUE, sep=",")
household_data2b_Nchelenge3 <- read.csv("ICEMRHouseholdMember_DATA.csv", header=TRUE, sep=",")

## went through each participant to examine which answers are most consistent, also looking at date of answer vs GPS
household_data2b_Nchelenge4 <- household_data2b_Nchelenge3[c(1,9,28,47,65,82,95,112,134,151,168,187,206,224,235,238,245,265,285,
                                                             293,305,314,328,349,354,360,361,380,399,404,409,414,419,421,440,459,478,
                                                             483,489,504,520,540,556,566,575,587,603,618,629,641,656,671,683,696,
                                                             708:722,723,735:739),]
household_data2b_Nchelenge5 <- household_data2b_Nchelenge4[,c(1,4:6,8)]


Nchelenge_temp <- Nchelenge[,c(11,6,7,16,34,28,32,33,41:45)]
Nchelenge_temp2 <- merge(Nchelenge5b2a2, Nchelenge_temp, by.x=c("partid","log_date", "log_time"), by.y=c("partid","date","time"), all.x=TRUE)
Nchelenge_temp2a <- Nchelenge_temp2[order(Nchelenge_temp2$partid, Nchelenge_temp2$point_num),]
Nchelenge_temp2b <- Nchelenge_temp2a[,c(1:57,60:64)]

Nchelenge_temp3 <- merge(Nchelenge_temp2b, household_data2b_Nchelenge5, by.x="partid", by.y="part_id", all.x=TRUE)
Nchelenge_temp3a <- Nchelenge_temp3[,c(1:54,63:66,58:62,55:57)]
##

#######
Nchelenge5_new <- Nchelenge5b2a2[,c(1:3,14,4,12,5:11,28,15:19,26,20:22,24,29,30,32,33,34:54,25)]
Nchelenge5_new2 <- Nchelenge5_new[order(Nchelenge5_new$partid, Nchelenge5_new$date_new, Nchelenge5_new$time_new_hour, Nchelenge5_new$time_new),]
Nchelenge_temp3b <- Nchelenge_temp3a[order(Nchelenge_temp3a$partid, Nchelenge_temp3a$date_new, Nchelenge_temp3a$time_new_hour, Nchelenge_temp3a$time_new),]

Nchelenge5_new3 <- cbind(Nchelenge5_new2,Nchelenge_temp3b[,c(55:66)])

##
############### CHANGE HOME CLUSTER FOR THOSE AT 2ND HOME ###############
Nchelenge5_new3$home_clust2_lat <- Nchelenge5_new3$hh_lat
Nchelenge5_new3$home_clust2_long <- Nchelenge5_new3$hh_long


# 21874 --> based on gps clusters, spend majority of time inland at farm location, whereas geopoint in lakeside
### so set household gps coordinates as lakeside (main house), but later set "home cluster" as farm house since closest hh dist to lake was 14km away (although most time at place 1km away)
### set home at place where most time spen in farm region
new_home_clust <- Nchelenge5_new3$clust[which(Nchelenge5_new3$partid=="21874")][which.max(Nchelenge5_new3$clust_trips_time_new[which(Nchelenge5_new3$partid=="21874")])]
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "21874" & Nchelenge5_new3$clust == new_home_clust)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "21874")] <- unique(Nchelenge5_new3$log_lat[which(Nchelenge5_new3$partid == "21874" & Nchelenge5_new3$home_clust2 == 1)])[1]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "21874")] <- unique(Nchelenge5_new3$log_long[which(Nchelenge5_new3$partid == "21874" & Nchelenge5_new3$home_clust2 == 1)])[1]


## 22117 --> have 2nd residence listed -->  -9.327  28.739 is lakeside (main home),  -9.363  28.873 is farm home
### set -9.327  28.739 as household gps here --> later set home cluster based on where spending time
## seems like home cluster is lakeside for 22118 and 22340 and farm house for 22117
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "22117" & Nchelenge5_new3$home_clust2 == 1)] <- 0 
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "22117" & Nchelenge5_new3$time_inclust > 10)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "22117")] <- unique(Nchelenge$hh_lat[which(Nchelenge$partid == "22117")])[1]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "22117")] <- unique(Nchelenge$hh_long[which(Nchelenge$partid == "22117")])[1]

## 22190 --> home listed is lakeside, but spending most time at other hh_lat/long listed, i.e., their farmhouse
# unique(Nchelenge$hh_lat[which(Nchelenge$partid == "22190")])
# unique(Nchelenge$hh_long[which(Nchelenge$partid == "22190")])
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "22190" & Nchelenge5_new3$home_clust2 == 1)] <- 0 
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "22190" & Nchelenge5_new3$time_inclust > 10)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "22190")] <- unique(Nchelenge$hh_lat[which(Nchelenge$partid == "22190")])[2]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "22190")] <- unique(Nchelenge$hh_long[which(Nchelenge$partid == "22190")])[2]

## 21897 --> have 2nd residence listed -->   -9.248521  28.80949 is lake (main home), -9.268372  28.84450 is farm home 
### so set household gps coordinates as lake (main house), but later set "home cluster" as farm house since most time spent there
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "21897" & Nchelenge5_new3$home_clust2 == 1)] <- 0 
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "21897" & Nchelenge5_new3$time_inclust > 10)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "21897")] <- unique(Nchelenge$hh_lat[which(Nchelenge$partid == "21897")])[1]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "21897")] <- unique(Nchelenge$hh_long[which(Nchelenge$partid == "21897")])[1]

## 21898 --> have 2nd residence listed -->   -9.248521  28.80949 is lake (main home), -9.268372  28.84450 is farm home 
### so set household gps coordinates as lake (main house), but later set "home cluster" as farm house since most time spent there
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "21898" & Nchelenge5_new3$home_clust2 == 1)] <- 0 
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "21898" & Nchelenge5_new3$time_inclust > 10)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "21898")] <- unique(Nchelenge$hh_lat[which(Nchelenge$partid == "21898")])[1]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "21898")] <- unique(Nchelenge$hh_long[which(Nchelenge$partid == "21898")])[1]

## 21863 --> have 2nd residence listed -->   -9.296998  28.7838 is farm (main home), -9.309166  28.7445 is lake home 
### so set household gps coordinates as farm (main house), but later set "home cluster" as lake house since most time spent there
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "21863" & Nchelenge5_new3$home_clust2 == 1)] <- 0 
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "21863" & Nchelenge5_new3$time_inclust > 10)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "21863")] <- unique(Nchelenge$hh_lat[which(Nchelenge$partid == "21863")])[2]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "21863")] <- unique(Nchelenge$hh_long[which(Nchelenge$partid == "21863")])[2]

## 22030 doesn't say has 2nd home, but has 2 hh_lat/long values listed in Nchelenge dataset, one lakeside
## 20 years old, male
unique(Nchelenge$hh_lat[which(Nchelenge$partid == "22030")])
unique(Nchelenge$hh_long[which(Nchelenge$partid == "22030")])
distHaversine(c(-9.290172, 28.75271), c(-9.309166, 28.7445)) ## 2km from where other hh coords listed, but still on lakeside
## spends 3 days at farmhouse, 5 hours at lake home, then 20 days at another home on lakeside
#### first few days, switches between farm home, lake home, then next couple weeks almost always staying at other place
#### list lakehome as main home cluster for gps for now
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "22030" & Nchelenge5_new3$home_clust2 == 1)] <- 0 
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "22030" & round(Nchelenge5_new3$log_lat, digits=3) == -9.309)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "22030")] <- unique(Nchelenge$hh_lat[which(Nchelenge$partid == "22030")])[2]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "22030")] <- unique(Nchelenge$hh_long[which(Nchelenge$partid == "22030")])[2]

## 22819 doesn't say has 2nd home, but has 2 hh_lat/long values listed in Nchelenge dataset, one lakeside
## 16 years old, female
unique(Nchelenge$hh_lat[which(Nchelenge$partid == "22819")])
unique(Nchelenge$hh_long[which(Nchelenge$partid == "22819")])
distHaversine(c(-9.355678,28.74494), c(-9.349862, 28.73764)) ## 1km from where other hh coords listed, but still on lakeside
## spends 1 days at farmhouse, 4 days at lake home, then 26 days at another home on lakeside
#### first few days at farm home, then a few days at lake home, then
## a couple weeks staying at other place, then a few days at lake home,
## then back to other place for 2 weeks
## other place gps not matched to another part's household gps
#### list lakehome as main home cluster for gps for now
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "22819" & Nchelenge5_new3$home_clust2 == 1)] <- 0 
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "22819" & round(Nchelenge5_new3$log_lat, digits=4) == -9.3556)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "22819")] <- unique(Nchelenge$hh_lat[which(Nchelenge$partid == "22819")])[2]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "22819")] <- unique(Nchelenge$hh_long[which(Nchelenge$partid == "22819")])[2]

## 22583 --> home listed is farmhouse, but spending time at lakeside
## spent 12 of 21 days lakeside. says has 2nd residence, but place is listed as Shondoni in southern zambia area, not lakeside
### didn't list trip
## 47 year old male; cluster 2 km from water
### —>—>spends a few hours at farm home, then takes about 30 minutes to go 5km back to lakeside and spends a few days at lakeside, and repeat.
### set home at place where most time spen in farm region
unique(Nchelenge$hh_lat[which(Nchelenge$partid == "22583")])
unique(Nchelenge$hh_long[which(Nchelenge$partid == "22583")])
new_home_clust <- Nchelenge5_new3$clust[which(Nchelenge5_new3$partid=="22583")][which.max(Nchelenge5_new3$clust_trips_time_new[which(Nchelenge5_new3$partid=="22583")])]
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "22583" & Nchelenge5_new3$home_clust2 == 1)] <- 0 
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "22583" & Nchelenge5_new3$clust == new_home_clust)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "22583")] <- unique(Nchelenge5_new3$log_lat[which(Nchelenge5_new3$partid == "22583" & Nchelenge5_new3$home_clust2 == 1)])[1]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "22583")] <- unique(Nchelenge5_new3$log_long[which(Nchelenge5_new3$partid == "22583" & Nchelenge5_new3$home_clust2 == 1)])[1]

## 21835
unique(Nchelenge$hh_lat[which(Nchelenge$partid == "21835")])
unique(Nchelenge$hh_long[which(Nchelenge$partid == "21835")])
Nchelenge5_new3[which(Nchelenge5_new3$partid == "21835" & Nchelenge5_new3$time_inclust > 10),][1,]
### 


## 22593
## 20 years old, female 
unique(Nchelenge$hh_lat[which(Nchelenge$partid == "22593")])
unique(Nchelenge$hh_long[which(Nchelenge$partid == "22593")])
Nchelenge5_new3[which(Nchelenge5_new3$partid == "22593" & Nchelenge5_new3$time_inclust > 10),][1,]
## spending 16 of 28 days 600 meters from listed home
## spends 8 of 28 days at home cluster (28%)
#### first few days at home, then a couple weeks staying at other place, 
## then a day at home, then back to other place for 1.5 weeks, then back home for a week
## other place gps not matched to another part's household gps
## other housemate is 19 year old female --> doesn't go to that house/gps area
#### -->--> for now, leave home same



a <- unique(Nchelenge5_new3$partid[which(Nchelenge5_new3$home_clust2 == 1)])
b <- unique(Nchelenge5_new3$partid)
no_home_clust <- b[which(!(b %in% a))]
### 21871 

## 21871
## all cluster >10km from home listed
unique(Nchelenge$hh_lat[which(Nchelenge$partid == "21871")])
unique(Nchelenge$hh_long[which(Nchelenge$partid == "21871")])
## only 1 listed
unique(household_data2b$part_id[which(household_data2b$hh_number == "204380")])
## 21874 housemate --> 
# 21874 --> based on gps clusters, spend majority of time inland at farm location, whereas geopoint in lakeside
### so set household gps coordinates as lakeside (main house), but later set "home cluster" as farm house since most time spent near there
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid=="21874")][1]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid=="21874")][1]

Nchelenge5_new3$log_lat[which(Nchelenge5_new3$partid == "21871")][which.max(Nchelenge5_new3$clust_trips_time_new[which(Nchelenge5_new3$partid=="21871")])]
Nchelenge5_new3$log_long[which(Nchelenge5_new3$partid == "21871")][which.max(Nchelenge5_new3$clust_trips_time_new[which(Nchelenge5_new3$partid=="21871")])]
## farm gps points are where 21871 spending most time near -->--> set farm house as home cluster for 21871 

new_home_clust <- Nchelenge5_new3$clust[which(Nchelenge5_new3$partid=="21871")][which.max(Nchelenge5_new3$clust_trips_time_new[which(Nchelenge5_new3$partid=="21871")])]

Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "21871" & Nchelenge5_new3$clust == new_home_clust)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "21871")] <- unique(Nchelenge5_new3$log_lat[which(Nchelenge5_new3$partid == "21871" & Nchelenge5_new3$home_clust2 == 1)])[1]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "21871")] <- unique(Nchelenge5_new3$log_long[which(Nchelenge5_new3$partid == "21871" & Nchelenge5_new3$home_clust2 == 1)])[1]

## 21890 --> have 2 hh-lat/hong values listed --> -9.3785 28.737 and -9.3882 28.736 --> both on lake, just slightly further down 
### so set household gps coordinates as place where most time spent since coordinates match 1st hh-lat/long
unique(Nchelenge$hh_lat[which(Nchelenge$partid == "21890")])
unique(Nchelenge$hh_long[which(Nchelenge$partid == "21890")])
Nchelenge5_new3$log_lat[which(Nchelenge5_new3$partid=="21890")][which.max(Nchelenge5_new3$clust_trips_time_new[which(Nchelenge5_new3$partid=="21890")])]
Nchelenge5_new3$log_long[which(Nchelenge5_new3$partid=="21890")][which.max(Nchelenge5_new3$clust_trips_time_new[which(Nchelenge5_new3$partid=="21890")])]

new_home_clust <- Nchelenge5_new3$clust[which(Nchelenge5_new3$partid=="21890")][which.max(Nchelenge5_new3$clust_trips_time_new[which(Nchelenge5_new3$partid=="21890")])]
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "21890" & Nchelenge5_new3$home_clust2 == 1)] <- 0 
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "21890" & Nchelenge5_new3$clust == new_home_clust)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "21890")] <- unique(Nchelenge5_new3$log_lat[which(Nchelenge5_new3$partid == "21890" & Nchelenge5_new3$home_clust2 == 1)])[1]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "21890")] <- unique(Nchelenge5_new3$log_long[which(Nchelenge5_new3$partid == "21890" & Nchelenge5_new3$home_clust2 == 1)])[1]


## 22581 --> no trips taken to place listed as home_clust, so removed --> ## place with most time is 19 meters form home loc (24 days spent)
## set new home_clust2 with this place
unique(Nchelenge$hh_lat[which(Nchelenge$partid == "22581")])
unique(Nchelenge$hh_long[which(Nchelenge$partid == "22581")])
Nchelenge5_new3$hhdist_m_haversine[which(Nchelenge5_new3$partid=="22581")][which.max(Nchelenge5_new3$clust_trips_time_new[which(Nchelenge5_new3$partid=="22581")])] # 19 meters

new_home_clust <- Nchelenge5_new3$clust[which(Nchelenge5_new3$partid=="22581")][which.max(Nchelenge5_new3$clust_trips_time_new[which(Nchelenge5_new3$partid=="22581")])]
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "21890" & Nchelenge5_new3$home_clust2 == 1)] <- 0 
Nchelenge5_new3$home_clust2[which(Nchelenge5_new3$partid == "22581" & Nchelenge5_new3$clust == new_home_clust)] <- 1
Nchelenge5_new3$home_clust2_lat[which(Nchelenge5_new3$partid == "22581")] <- unique(Nchelenge5_new3$log_lat[which(Nchelenge5_new3$partid == "22581" & Nchelenge5_new3$home_clust2 == 1)])[1]
Nchelenge5_new3$home_clust2_long[which(Nchelenge5_new3$partid == "22581")] <- unique(Nchelenge5_new3$log_long[which(Nchelenge5_new3$partid == "22581" & Nchelenge5_new3$home_clust2 == 1)])[1]
### 

##
#### redefine hh_dist for new home_clust2 ####
Nchelenge5_new3$hhdist_m_haversine2<- distHaversine(Nchelenge5_new3[,11:12], Nchelenge5_new3[,63:64])

############# fix those with points in lake ##############
lake_bound <- st_read("Lake_Mweru2.shp")
lake_bound2 <- as_Spatial(lake_bound)

Nchelenge5_new3_one_each <- Nchelenge5_new3[!duplicated(Nchelenge5_new3[,c(1,21)]),]

point <- data.frame(lon=Nchelenge5_new3_one_each$log_long, lat=Nchelenge5_new3_one_each$log_lat)
sp2   <- SpatialPoints(point,proj4string=CRS(proj4string(lake_bound2)))
Nchelenge5_new3_one_each$in_lake <- as.factor(as.integer(gContains(lake_bound2,sp2, byid=TRUE)))
#####
Nchelenge5_new3_one_each$lake_dist[which(Nchelenge5_new3_one_each$in_lake == 1 & Nchelenge5_new3_one_each$partid == "22816")]
##1480 meters from edge of lake --> other points near to lake edge --> find point closest to this cluster
dists <- distHaversine(Nchelenge5_new3_one_each[which(Nchelenge5_new3_one_each$partid == "22816"),11:12], 
                       Nchelenge5_new3_one_each[which(Nchelenge5_new3_one_each$partid == "22816" & Nchelenge5_new3_one_each$in_lake == 1),11:12])
dist1 <- order(dists)[2]
ind <- which(Nchelenge5_new3_one_each$partid == "22816")[dist1] #684
Nchelenge5_new3_one_each$lake_dist[ind] #157 meters from lake

ind2 <- which(Nchelenge5_new3_one_each$partid == "22816" & Nchelenge5_new3_one_each$in_lake == 1) 
clust_keep <- Nchelenge5_new3_one_each$clust[ind]
clust_rem <- Nchelenge5_new3_one_each$clust[ind2]

Nchelenge5_new3$clust[which(Nchelenge5_new3$partid=="22816" & Nchelenge5_new3$clust == clust_rem)] <- clust_keep
Nchelenge5_new3$points_inclust[which(Nchelenge5_new3$partid=="22816" & Nchelenge5_new3$clust == clust_keep)] <- Nchelenge5_new3_one_each$points_inclust[ind] + Nchelenge5_new3_one_each$points_inclust[ind2]
Nchelenge5_new3$time_inclust[which(Nchelenge5_new3$partid=="22816" & Nchelenge5_new3$clust == clust_keep)] <- Nchelenge5_new3_one_each$time_inclust[ind] + Nchelenge5_new3_one_each$time_inclust[ind2]
Nchelenge5_new3$time_inclust_min[which(Nchelenge5_new3$partid=="22816" & Nchelenge5_new3$clust == clust_keep)] <- Nchelenge5_new3_one_each$time_inclust_min[ind] + Nchelenge5_new3_one_each$time_inclust_min[ind2]

Nchelenge5_new3$clust_trips_time_new[which(Nchelenge5_new3$partid=="22816" & Nchelenge5_new3$clust==clust_keep)] <- Nchelenge5_new3_one_each$clust_trips_time_new[ind] + Nchelenge5_new3_one_each$clust_trips_time_new[ind2]
Nchelenge5_new3$clust_trips_time_new_morning[which(Nchelenge5_new3$partid=="22816" & Nchelenge5_new3$clust==clust_keep)] <- Nchelenge5_new3_one_each$clust_trips_time_new_morning[ind] + Nchelenge5_new3_one_each$clust_trips_time_new_morning[ind2]
Nchelenge5_new3$clust_trips_time_new_midday[which(Nchelenge5_new3$partid=="22816" & Nchelenge5_new3$clust==clust_keep)] <- Nchelenge5_new3_one_each$clust_trips_time_new_midday[ind] + Nchelenge5_new3_one_each$clust_trips_time_new_midday[ind2]
Nchelenge5_new3$clust_trips_time_new_evening[which(Nchelenge5_new3$partid=="22816" & Nchelenge5_new3$clust==clust_keep)] <- Nchelenge5_new3_one_each$clust_trips_time_new_evening[ind] + Nchelenge5_new3_one_each$clust_trips_time_new_evening[ind2]
Nchelenge5_new3$clust_trips_time_new_night[which(Nchelenge5_new3$partid=="22816" & Nchelenge5_new3$clust==clust_keep)] <- Nchelenge5_new3_one_each$clust_trips_time_new_night[ind] + Nchelenge5_new3_one_each$clust_trips_time_new_night[ind2]

inds <- which(Nchelenge5_new3$partid=="22816" & Nchelenge5_new3$clust==clust_keep)
counts <- cumsum(c(1,diff(inds)>6))
counts[which(counts %in% c(which(summary(as.factor(counts), maxsum=500) < 7)))]<- NA
counts <- as.factor(counts)
levels(counts) <- c(1:nlevels(counts))

Nchelenge5_new3$trip_count_new[which(Nchelenge5_new3$partid=="22816" & Nchelenge5_new3$clust==clust_keep)] <- as.numeric(as.character(counts))
Nchelenge5_new3_temp <- Nchelenge5_new3[inds,]
xx <- aggregate(time_atpoint ~ trip_count_new, Nchelenge5_new3_temp, sum)
for(kk in 1:nrow(xx)){
  Nchelenge5_new3$trip_time_new[which(Nchelenge5_new3$partid=="22816" &Nchelenge5_new3$clust == clust_keep & Nchelenge5_new3$trip_count_new == kk)] <- xx$time_atpoint[kk]
}
Nchelenge5_new3$trip_time_new[which(is.na(Nchelenge5_new3$trip_count_new))] <- NA
Nchelenge5_new3$trip_time_new_min[inds] <- Nchelenge5_new3$trip_time_new[inds]*24*60

Nchelenge5_new3b <- Nchelenge5_new3[order(Nchelenge5_new3$partid, Nchelenge5_new3$date_new, Nchelenge5_new3$time_new_hour, Nchelenge5_new3$time_new),]

trips <- levels(as.factor(Nchelenge5_new3b$trip_count_new[which(Nchelenge5_new3b$partid == "22816" & Nchelenge5_new3b$clust == clust_keep)]))
for(kk in 1:length(trips)){
  trip <- trips[kk]
  temp2 <- Nchelenge5_new3b[which(Nchelenge5_new3b$partid == "22816" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == trip),]
  time_at_points_middle <- diff(c(as.POSIXct(paste(temp2$date_new, temp2$time_new), format="%Y-%m-%d %H:%M:%S", tz="Africa/Lusaka")),2)
  if(units(time_at_points_middle) == "hours"){
    time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/24
  }
  if(units(time_at_points_middle) == "mins"){
    time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/60/24
  }
  if(units(time_at_points_middle) == "secs"){
    time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/60/60/24
  }
  time_at_points <- c(as.numeric(difftime(as.POSIXct(paste(temp2$date_new[2], temp2$time_new[2])), as.POSIXct(paste(temp2$date_new[1], temp2$time_new[1])), unit="min"))/2/60/24, 
                      time_at_points_middle2,
                      as.numeric(difftime(as.POSIXct(paste(temp2$date_new[c(nrow(temp2))], temp2$time_new[c(nrow(temp2))])), 
                                          as.POSIXct(paste(temp2$date_new[c(nrow(temp2)-1)], temp2$time_new[c(nrow(temp2)-1)])), 
                                          unit="min"))/2/60/24)
  Nchelenge5_new3b$time_atpoint_new[which(Nchelenge5_new3b$partid == "22816" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == trip)] <- time_at_points
}
temp3 <- Nchelenge5_new3b[which(Nchelenge5_new3b$partid == "22816" & Nchelenge5_new3b$clust == clust_keep),]
xx <- aggregate(time_atpoint_new ~ trip_count_new, temp3, sum)
for(kk in 1:nrow(xx)){
  Nchelenge5_new3b$trip_time_new2[which(Nchelenge5_new3b$partid == "22816" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == kk)] <- xx$time_atpoint[kk]
}
xx2 <- aggregate(time_atpoint_new ~ trip_count_new + time_of_day, temp3, sum)
xx2_morning <- xx2[which(xx2$time_of_day == "morning"),]
xx2_midday <- xx2[which(xx2$time_of_day == "midday"),]
xx2_evening <- xx2[which(xx2$time_of_day == "evening"),]
xx2_night <- xx2[which(xx2$time_of_day == "night"),]
for(kk in 1:nrow(xx2_morning)){
  Nchelenge5_new3b$trip_time_new2_morning[which(Nchelenge5_new3b$partid == "22816" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == xx2_morning$trip_count_new[kk])] <- xx2_morning$time_atpoint_new[kk]
}
for(kk in 1:nrow(xx2_midday)){
  Nchelenge5_new3b$trip_time_new2_midday[which(Nchelenge5_new3b$partid == "22816" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == xx2_midday$trip_count_new[kk])] <- xx2_midday$time_atpoint_new[kk]
}
for(kk in 1:nrow(xx2_evening)){
  Nchelenge5_new3b$trip_time_new2_evening[which(Nchelenge5_new3b$partid == "22816" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == xx2_evening$trip_count_new[kk])] <- xx2_evening$time_atpoint_new[kk]
}
for(kk in 1:nrow(xx2_night)){
print(xx2_night$trip_count_new[kk])
  }
##
######
sum(Nchelenge5_new3_one_each$time_inclust[which(Nchelenge5_new3_one_each$partid == "22029" & Nchelenge5_new3_one_each$in_lake == 0)])
sum(Nchelenge5_new3_one_each$time_inclust[which(Nchelenge5_new3_one_each$partid == "22029" & Nchelenge5_new3_one_each$in_lake == 1)])
## spends 27 days out of lake and 3 days in lake --> other points near to lake edge --> find point closest to this cluster
## closest one
which.min(Nchelenge5_new3_one_each$lake_dist[which(Nchelenge5_new3_one_each$partid == "22029" & Nchelenge5_new3_one_each$in_lake == 1)])
ind <- which(Nchelenge5_new3_one_each$partid == "22029" & Nchelenge5_new3_one_each$in_lake == 1)[which.min(Nchelenge5_new3_one_each$lake_dist[which(Nchelenge5_new3_one_each$partid == "22029" & Nchelenge5_new3_one_each$in_lake == 1)])]
dists <- distHaversine(Nchelenge5_new3_one_each[ind,11:12], 
                       Nchelenge5_new3_one_each[which(Nchelenge5_new3_one_each$partid == "22029" & Nchelenge5_new3_one_each$in_lake == 0),11:12])
dist1 <- order(dists)[1]
ind2 <- which(Nchelenge5_new3_one_each$partid == "22029" & Nchelenge5_new3_one_each$in_lake == 0)[dist1]
Nchelenge5_new3_one_each$lake_dist[ind2] #110 meters from lake
ind3 <- which(Nchelenge5_new3_one_each$partid == "22029" & Nchelenge5_new3_one_each$in_lake == 1)

clust_keep <- Nchelenge5_new3_one_each$clust[ind2]
clusts_rem <- Nchelenge5_new3_one_each$clust[which(Nchelenge5_new3_one_each$partid == "22029" & Nchelenge5_new3_one_each$in_lake == 1)]

### add time at lake points to this cluster
Nchelenge5_new3b$clust[which(Nchelenge5_new3b$partid=="22029" & Nchelenge5_new3b$clust %in% clusts_rem)] <- clust_keep
Nchelenge5_new3b$points_inclust[which(Nchelenge5_new3b$partid=="22029" & Nchelenge5_new3b$clust == clust_keep)] <- sum(Nchelenge5_new3_one_each$points_inclust[ind2], Nchelenge5_new3_one_each$points_inclust[ind3])
Nchelenge5_new3b$time_inclust[which(Nchelenge5_new3b$partid=="22029" & Nchelenge5_new3b$clust == clust_keep)] <- sum(Nchelenge5_new3_one_each$time_inclust[ind2], Nchelenge5_new3_one_each$time_inclust[ind3])
Nchelenge5_new3b$time_inclust_min[which(Nchelenge5_new3b$partid=="22029" & Nchelenge5_new3b$clust == clust_keep)] <-sum(Nchelenge5_new3_one_each$time_inclust_min[ind2], Nchelenge5_new3_one_each$time_inclust_min[ind3])

Nchelenge5_new3b$clust_trips_time_new[which(Nchelenge5_new3b$partid=="22029" & Nchelenge5_new3b$clust==clust_keep)] <- sum(Nchelenge5_new3_one_each$clust_trips_time_new[ind2], Nchelenge5_new3_one_each$clust_trips_time_new[ind3])
Nchelenge5_new3b$clust_trips_time_new_morning[which(Nchelenge5_new3b$partid=="22029" & Nchelenge5_new3b$clust==clust_keep)] <- sum(Nchelenge5_new3_one_each$clust_trips_time_new_morning[ind2], Nchelenge5_new3_one_each$clust_trips_time_new_morning[ind3], na.rm=TRUE)
Nchelenge5_new3b$clust_trips_time_new_midday[which(Nchelenge5_new3b$partid=="22029" & Nchelenge5_new3b$clust==clust_keep)] <- sum(Nchelenge5_new3_one_each$clust_trips_time_new_midday[ind2], Nchelenge5_new3_one_each$clust_trips_time_new_midday[ind3], na.rm=TRUE)
Nchelenge5_new3b$clust_trips_time_new_evening[which(Nchelenge5_new3b$partid=="22029" & Nchelenge5_new3b$clust==clust_keep)] <- sum(Nchelenge5_new3_one_each$clust_trips_time_new_evening[ind2], Nchelenge5_new3_one_each$clust_trips_time_new_evening[ind3], na.rm=TRUE)
Nchelenge5_new3b$clust_trips_time_new_night[which(Nchelenge5_new3b$partid=="22029" & Nchelenge5_new3b$clust==clust_keep)] <- sum(Nchelenge5_new3_one_each$clust_trips_time_new_night[ind2], Nchelenge5_new3_one_each$clust_trips_time_new_night[ind3], na.rm=TRUE)

inds <- which(Nchelenge5_new3b$partid=="22029" & Nchelenge5_new3b$clust==clust_keep)
counts <- cumsum(c(1,diff(inds)>6))
counts[which(counts %in% c(which(summary(as.factor(counts), maxsum=500) < 7)))]<- NA
counts <- as.factor(counts)
levels(counts) <- c(1:nlevels(counts))

Nchelenge5_new3b$trip_count_new[which(Nchelenge5_new3b$partid=="22029" & Nchelenge5_new3b$clust==clust_keep)] <- as.numeric(as.character(counts))
Nchelenge5_new3b_temp <- Nchelenge5_new3b[inds,]
xx <- aggregate(time_atpoint ~ trip_count_new, Nchelenge5_new3b_temp, sum)
for(kk in 1:nrow(xx)){
  Nchelenge5_new3b$trip_time_new[which(Nchelenge5_new3b$partid=="22029" &Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == kk)] <- xx$time_atpoint[kk]
}
Nchelenge5_new3b$trip_time_new[which(is.na(Nchelenge5_new3b$trip_count_new))] <- NA
Nchelenge5_new3b$trip_time_new_min[inds] <- Nchelenge5_new3b$trip_time_new[inds]*24*60

Nchelenge5_new3b <- Nchelenge5_new3b[order(Nchelenge5_new3b$partid, Nchelenge5_new3b$date_new, Nchelenge5_new3b$time_new_hour, Nchelenge5_new3b$time_new),]

trips <- levels(as.factor(Nchelenge5_new3b$trip_count_new[which(Nchelenge5_new3b$partid == "22029" & Nchelenge5_new3b$clust == clust_keep)]))
for(kk in 1:length(trips)){
  trip <- trips[kk]
  temp2 <- Nchelenge5_new3b[which(Nchelenge5_new3b$partid == "22029" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == trip),]
  time_at_points_middle <- diff(c(as.POSIXct(paste(temp2$date_new, temp2$time_new), format="%Y-%m-%d %H:%M:%S", tz="Africa/Lusaka")),2)
  if(units(time_at_points_middle) == "hours"){
    time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/24
  }
  if(units(time_at_points_middle) == "mins"){
    time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/60/24
  }
  if(units(time_at_points_middle) == "secs"){
    time_at_points_middle2 <- as.numeric(time_at_points_middle)/2/60/60/24
  }
  time_at_points <- c(as.numeric(difftime(as.POSIXct(paste(temp2$date_new[2], temp2$time_new[2])), as.POSIXct(paste(temp2$date_new[1], temp2$time_new[1])), unit="min"))/2/60/24, 
                      time_at_points_middle2,
                      as.numeric(difftime(as.POSIXct(paste(temp2$date_new[c(nrow(temp2))], temp2$time_new[c(nrow(temp2))])), 
                                          as.POSIXct(paste(temp2$date_new[c(nrow(temp2)-1)], temp2$time_new[c(nrow(temp2)-1)])), 
                                          unit="min"))/2/60/24)
  Nchelenge5_new3b$time_atpoint_new[which(Nchelenge5_new3b$partid == "22029" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == trip)] <- time_at_points
}
temp3 <- Nchelenge5_new3b[which(Nchelenge5_new3b$partid == "22029" & Nchelenge5_new3b$clust == clust_keep),]
xx <- aggregate(time_atpoint_new ~ trip_count_new, temp3, sum)
for(kk in 1:nrow(xx)){
  Nchelenge5_new3b$trip_time_new2[which(Nchelenge5_new3b$partid == "22029" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == kk)] <- xx$time_atpoint[kk]
}
xx2 <- aggregate(time_atpoint_new ~ trip_count_new + time_of_day, temp3, sum)
xx2_morning <- xx2[which(xx2$time_of_day == "morning"),]
xx2_midday <- xx2[which(xx2$time_of_day == "midday"),]
xx2_evening <- xx2[which(xx2$time_of_day == "evening"),]
xx2_night <- xx2[which(xx2$time_of_day == "night"),]
for(kk in 1:nrow(xx2_morning)){
  Nchelenge5_new3b$trip_time_new2_morning[which(Nchelenge5_new3b$partid == "22029" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == xx2_morning$trip_count_new[kk])] <- xx2_morning$time_atpoint_new[kk]
}
for(kk in 1:nrow(xx2_midday)){
  Nchelenge5_new3b$trip_time_new2_midday[which(Nchelenge5_new3b$partid == "22029" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == xx2_midday$trip_count_new[kk])] <- xx2_midday$time_atpoint_new[kk]
}
for(kk in 1:nrow(xx2_evening)){
  Nchelenge5_new3b$trip_time_new2_evening[which(Nchelenge5_new3b$partid == "22029" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == xx2_evening$trip_count_new[kk])] <- xx2_evening$time_atpoint_new[kk]
}
for(kk in 1:nrow(xx2_night)){
  Nchelenge5_new3b$trip_time_new2_night[which(Nchelenge5_new3b$partid == "22029" & Nchelenge5_new3b$clust == clust_keep & Nchelenge5_new3b$trip_count_new == xx2_night$trip_count_new[kk])] <- xx2_night$time_atpoint_new[kk]
}

Nchelenge5_new3b_one_each <- Nchelenge5_new3b[!duplicated(Nchelenge5_new3b[,c(1,21)]),]
point <- data.frame(lon=Nchelenge5_new3b_one_each$log_long, lat=Nchelenge5_new3b_one_each$log_lat)
sp2   <- SpatialPoints(point,proj4string=CRS(proj4string(lake_bound2)))
Nchelenge5_new3b_one_each$in_lake <- as.factor(as.integer(gContains(lake_bound2,sp2, byid=TRUE)))

summary(Nchelenge5_new3b_one_each$in_lake) ## 1 left --> guy 65 meters in 
##

############ 
############ 
Nchelenge5_new3b2 <- Nchelenge5_new3b[order(Nchelenge5_new3b$partid, Nchelenge5_new3b$date_new, Nchelenge5_new3b$time_new_hour, Nchelenge5_new3b$time_new),]

### take out those with <15 minutes in cluster or <15 minute trips that may be here now because of clusters that were >1000 meter radius
Nchelenge5_new3b22 <- Nchelenge5_new3b2[which(Nchelenge5_new3b2$clust_trips_time_new > 0.01),]
Nchelenge5_new3b22b <- Nchelenge5_new3b22[which(Nchelenge5_new3b22$trip_time_new2 > 0.01),]


Nchelenge5_clusts_shorter <- Nchelenge5_new3b22b[order(Nchelenge5_new3b22b$partid, Nchelenge5_new3b22b$clust),]
Nchelenge5_clusts_shorter2 <- Nchelenge5_clusts_shorter[!duplicated(Nchelenge5_clusts_shorter[,c(1,21)]),c(1:6,15,16,11,12,18:24,28,29:34,45:49,40:44,50,51:65)]

Nchelenge5_clust_trips_shorter <- Nchelenge5_new3b22b[order(Nchelenge5_new3b22b$partid, Nchelenge5_new3b22b$clust, Nchelenge5_new3b22b$trip_count_new),]
Nchelenge5_clust_trips_shorter2 <- Nchelenge5_clust_trips_shorter[!duplicated(Nchelenge5_clust_trips_shorter[,c(1,21,25)]),c(1:6,15,16,11,12,18:22,23,25:27,28,29:34,45:49,40:44,35:39,50,51:65)]

##
#####
setwd()
write.csv(Nchelenge4, "Nchelenge4.csv", row.names = FALSE)
write.csv(Nchelenge4b, "Nchelenge4b.csv", row.names = FALSE)
write.csv(Nchelenge5, "Nchelenge5.csv", row.names = FALSE)
write.csv(Nchelenge5_new3b22b, "Nchelenge5_new.csv", row.names = FALSE)

write.csv(Nchelenge5_clusts_shorter2, "Nchelenge5_clusts_shorter2.csv", row.names = FALSE)
write.csv(Nchelenge5_clust_trips_shorter2, "Nchelenge5_clust_trips_shorter2.csv", row.names = FALSE)
##
############ 