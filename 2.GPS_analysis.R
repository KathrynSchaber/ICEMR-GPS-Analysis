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
library(sf)
library(plotly)
library(ggpubr)
library(fitdistrplus)
library(ggplot2)
library(lme4)
library(bbmle)
library(emmeans)
library(ggridges)
library(gamlss)
library(corrplot)
library(PerformanceAnalytics)
library(Hmisc)
library(colorBlindness)
##### read in datasets ####
setwd()
Mutasa5 <- read.csv("Mutasa5.csv", header=TRUE, sep=",") # all points
Choma5 <- read.csv("Choma5.csv", header=TRUE, sep=",") # all points
Nchelenge5 <- read.csv("Nchelenge5.csv", header=TRUE, sep=",") # all points

Mutasa5_new <- read.csv("Mutasa5_new.csv", header=TRUE, sep=",") # only points with keep=1 (no outliers, only top 15 clusters by time or points spent, no clusters with <5 points or <30 minutes spent)
Mutasa5_clusts_shorter2 <- read.csv("Mutasa5_clusts_shorter2.csv", header=TRUE, sep=",") # one line per cluster
Mutasa5_clust_trips_shorter2 <- read.csv("Mutasa5_clust_trips_shorter2.csv", header=TRUE, sep=",") # one line per trip to cluster

Choma5_new <- read.csv("Choma5_new.csv", header=TRUE, sep=",") # only points with keep=1 (no outliers, only top 15 clusters by time or points spent, no clusters with <5 points or <30 minutes spent)
Choma5_clusts_shorter2 <- read.csv("Choma5_clusts_shorter2.csv", header=TRUE, sep=",") # one line per cluster
Choma5_clust_trips_shorter2 <- read.csv("Choma5_clust_trips_shorter2.csv", header=TRUE, sep=",") # one line per trip to cluster

Nchelenge5_new <- read.csv("Nchelenge5_new.csv", header=TRUE, sep=",") # only points with keep=1 (no outliers, only top 15 clusters by time or points spent, no clusters with <5 points or <30 minutes spent)
Nchelenge5_clusts_shorter2 <- read.csv("Nchelenge5_clusts_shorter2.csv", header=TRUE, sep=",") # one line per cluster
Nchelenge5_clust_trips_shorter2 <- read.csv("Nchelenge5_clust_trips_shorter2.csv", header=TRUE, sep=",") # one line per trip to cluster
##
###############
###############
##### CHOMA #####
setwd("/Users/kschaber/Desktop/ICEMR_GPSLoggers/Data")
zambia <- st_read("Zambia.shp")
zambia2 <- as_Spatial(zambia)

point <- data.frame(lon=Choma5_clusts_shorter2$log_long, lat=Choma5_clusts_shorter2$log_lat)
sp2   <- SpatialPoints(point,proj4string=CRS(proj4string(zambia2)))
Choma5_clusts_shorter2$in_zambia <- as.factor(as.integer(gContains(zambia2,sp2, byid=TRUE)))
setwd("/Users/kschaber/Desktop/ICEMR_GPSLoggers/Code/HDBSCAN/Outputs/New/3.17")

summary(Choma5_clusts_shorter2$in_zambia) ## no one leaves country

#### Get values for biting time period ####
Choma5_new2 <- Choma5_new
Choma5_new2$biting_time <- factor(0, levels = c(0,1))
Choma5_new2$biting_time[which(Choma5_new2$time_new_hour %in% c(21:23,0:5))] <- 1

Choma5_new2$trip_time_new2_biting_time <- as.numeric(NA)
Choma5_new2$clust_trips_time_new_biting_time <- as.numeric(NA)

parts <- levels(as.factor(Choma5_new2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Choma5_new2$clust[which(Choma5_new2$partid == part)]))
  for(jj in 1:length(clusts)){
    clust <- clusts[jj]
    trips <- levels(as.factor(Choma5_new2$trip_count_new[which(Choma5_new2$partid == part & Choma5_new2$clust == clust)]))
    temp3 <- Choma5_new2[which(Choma5_new2$partid == part & Choma5_new2$clust == clust & Choma5_new2$biting_time == 1),]
    if(nrow(temp3) > 0){
      xx2 <- aggregate(time_atpoint_new ~ trip_count_new, temp3, sum)
      for(kk in 1:nrow(xx2)){
        Choma5_new2$trip_time_new2_biting_time[which(Choma5_new2$partid == part & 
                                                        Choma5_new2$clust == clust &
                                                        Choma5_new2$trip_count_new == xx2$trip_count_new[kk])] <- xx2$time_atpoint_new[kk]
      }
    }
  }
  temp4 <- Choma5_new2[which(Choma5_new2$partid == part),]
  temp4b <- temp4[!duplicated(temp4[c(1,22,26)]),]
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_biting_time))),])>0){
    biting <- aggregate(trip_time_new2_biting_time ~ clust, temp4b, sum)
    if(nrow(biting) > 0){
      for(kk in 1:nrow(biting)){
        Choma5_new2$clust_trips_time_new_biting_time[which(Choma5_new2$partid == part & Choma5_new2$clust == biting$clust[kk] & !(is.na(Choma5_new2$trip_count_new)))] <- biting$trip_time_new2_biting_time[kk]
      }
    }
  }
}
Choma5_new2$part_trips_time_new_biting_time <- as.numeric(0)
temp <- Choma5_new2[!duplicated(Choma5_new2[c(1,22)]),]
tot <- aggregate(clust_trips_time_new_biting_time ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Choma5_new2$part_trips_time_new_biting_time[which(Choma5_new2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_biting_time[ii]
}

prop.table(table(Choma5_new2$home_clust2, Choma5_new2$biting_time, deparse.level = 2), margin=2)*100
#                       Choma5_new2$biting_time
# Choma5_new2$home_clust2           0           1
# 0                       58.64968174 50.28852737
# 1                       41.35031826 49.71147263
prop.table(table(Choma5_new2$home_clust2, Choma5_new2$biting_time, deparse.level = 2), margin=1)*100
#                       Choma5_new2$biting_time
# Choma5_new2$home_clust2           0           1
# 0                       82.25502151 17.74497849
# 1                       76.77698411 23.22301589

#### 50% of observations taken during biting times were in the home cluster
#### 23% of obervations in the home cluster were taken during biting times

Choma5_clusts_shorter_temp <- Choma5_new2[order(Choma5_new2$partid, Choma5_new2$clust),]
Choma5_clusts_shorter_temp2 <- Choma5_clusts_shorter_temp[!duplicated(Choma5_clusts_shorter_temp[,c(1,22)]),c(1:7,16,17,12,13,19:25,29,30:35,46:50,41:45,51,55,54)]

Choma5_clust_trips_shorter_temp <- Choma5_new2[order(Choma5_new2$partid, Choma5_new2$clust, Choma5_new2$trip_count),]
Choma5_clust_trips_shorter_temp2 <- Choma5_clust_trips_shorter_temp[!duplicated(Choma5_clust_trips_shorter_temp[,c(1,22,26)]),c(1:7,16,17,12,13,19:22,24,26:28,29,30:35,46:50,41:45,36:40,51,55,54,53)]

Choma5_clusts_shorter2$part_trips_time_new_biting_time <- Choma5_clusts_shorter_temp2$part_trips_time_new_biting_time
Choma5_clusts_shorter2$clust_trips_time_new_biting_time <- Choma5_clusts_shorter_temp2$clust_trips_time_new_biting_time

Choma5_clust_trips_shorter2$part_trips_time_new_biting_time <- Choma5_clust_trips_shorter_temp2$part_trips_time_new_biting_time
Choma5_clust_trips_shorter2$clust_trips_time_new_biting_time <- Choma5_clust_trips_shorter_temp2$clust_trips_time_new_biting_time
Choma5_clust_trips_shorter2$trip_time_new2_biting_time <- Choma5_clust_trips_shorter_temp2$trip_time_new2_biting_time
##
### number of locations (without home) ####
Choma5_clusts_shorter2_no_home <- Choma5_clusts_shorter2[which(Choma5_clusts_shorter2$home_clust2 == 0),]
loc_counts <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)$n
only_home <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home$partid))
loc_counts2 <- c(loc_counts,rep(0, only_home))
favstats(loc_counts2)
# min Q1 median Q3 max     mean       sd  n missing
#    0 7.25   12 17  22 11.90323 5.730876 62       0
hist(loc_counts2, xlim=c(0,26))
#
### distance of locations (without home) ####
### TOTAL ###
Choma5_clusts_shorter2_no_home$hhdist_km_haversine <- Choma5_clusts_shorter2_no_home$hhdist_m_haversine/1000
favstats(Choma5_clusts_shorter2_no_home$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
# 0.03054886 1.102931 7.930858 15.90348 112.2062 10.82367 13.58644 738       0

### PER PART ###
medians_km <- (favstats(Choma5_clusts_shorter2_no_home$hhdist_m_haversine~Choma5_clusts_shorter2_no_home$partid)$median)/1000
favstats(medians_km) # km
#          min        Q1   median       Q3      max     mean       sd  n missing
#  0.2448758 0.9599718 5.602667 14.34616 32.16992 8.109393 7.934682 61       0

### PLOT BOTH ###
x <- seq(0,500, length.out=500)
df <- with(Choma5_clusts_shorter2_no_home, data.frame(x = x, y1 = dnorm(x, mean(hhdist_km_haversine), sd(hhdist_km_haversine)),
                                                      y2 = dnorm(x, mean(medians_km), sd(medians_km))))

dists_km2 <- as.data.frame(cbind(c(Choma5_clusts_shorter2_no_home$hhdist_km_haversine, medians_km),
                                 c(rep("overall",nrow(Choma5_clusts_shorter2_no_home)),
                                   rep("part_median",length(medians_km)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,120))

## just medians
ggplot(data=dists_km2[which(dists_km2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 5, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,35))
#
### number of trips (without home) ####
### TOTAL ###
Choma5_clust_trips_shorter2_no_home <- Choma5_clust_trips_shorter2[which(Choma5_clust_trips_shorter2$home_clust2 == 0),]
Choma5_clust_trips_shorter2_no_home2 <- Choma5_clust_trips_shorter2_no_home[order(Choma5_clust_trips_shorter2_no_home$partid,
                                                                                  Choma5_clust_trips_shorter2_no_home$clust,
                                                                                  -Choma5_clust_trips_shorter2_no_home$trip_count_new),]
Choma5_clust_trips_shorter2_no_home3 <- Choma5_clust_trips_shorter2_no_home2[!duplicated(Choma5_clust_trips_shorter2_no_home2[,c(1,15)]),]

favstats(Choma5_clust_trips_shorter2_no_home3$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      2  3  58 3.063686 5.280209 738       0

### PER PART ###
trips_part_median <- favstats(Choma5_clust_trips_shorter2_no_home3$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3$partid)$median
favstats(trips_part_median)
#   min Q1 median Q3 max     mean       sd  n missing
#    1  1    1  2   4 1.540984 0.6664276 61       0

### PLOT BOTH ###
x <- seq(0,500, length.out=1000)
df <- with(Choma5_clust_trips_shorter2_no_home3, data.frame(x = x, y1 = dnorm(x, mean(trip_count_new), sd(trip_count_new)),
                                                            y2 = dnorm(x, mean(trips_part_median), sd(trips_part_median))))

trips_part2 <- as.data.frame(cbind(c(Choma5_clust_trips_shorter2_no_home3$trip_count_new, trips_part_median),
                                   c(rep("overall",nrow(Choma5_clust_trips_shorter2_no_home3)),
                                     rep("part_median",length(trips_part_median)))))
colnames(trips_part2) <- c("values","type")
trips_part2$values <- as.numeric(as.character(trips_part2$values))
trips_part2$type <- factor(trips_part2$type)

ggplot(data=trips_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,60))

## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1.5) +
  coord_cartesian(xlim=c(0,10))+ 
  scale_x_continuous(breaks=c(seq(0,10,2)))
#
### time per trip ####
### TOTAL ###
Choma5_clust_trips_shorter2_no_home$trip_time_new_hour <- Choma5_clust_trips_shorter2_no_home$trip_time_new2*24
favstats(Choma5_clust_trips_shorter2_no_home$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2402778 0.6800694   1.44 4.099028 472.56 9.028822 32.20347 2236       0

### PER PART ###
trip_time_part_median <- favstats(Choma5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home$partid)$median
favstats(trip_time_part_median)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.5879167 1.173333   1.44 1.748472 5.650556 1.636708 0.9599658 61       0

### PLOT BOTH ###
x <- seq(0,500, length.out=500)
df <- with(Choma5_clust_trips_shorter2_no_home, data.frame(x = x, y1 = dnorm(x, mean(trip_time_new_hour), sd(trip_time_new_hour)),
                                                           y2 = dnorm(x, mean(trip_time_part_median), sd(trip_time_part_median))))

trip_times_part2 <- as.data.frame(cbind(c(Choma5_clust_trips_shorter2_no_home$trip_time_new_hour, trip_time_part_median),
                                        c(rep("overall",nrow(Choma5_clust_trips_shorter2_no_home)),
                                          rep("part_median",length(trip_time_part_median)))))
colnames(trip_times_part2) <- c("values","type")
trip_times_part2$values <- as.numeric(as.character(trip_times_part2$values))
trip_times_part2$type <- factor(trip_times_part2$type)

ggplot(data=trip_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,500))


## just medians
ggplot(data=trip_times_part2[which(trip_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 0.5, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,6))
#
#
### Home Loc ####
Choma5_clusts_shorter2_home <- Choma5_clusts_shorter2[which(Choma5_clusts_shorter2$home_clust2 == 1),]
no_home <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_home$partid))# 4
## 3 people with no homes 
# # "10692b" "10753" "15927"
### number of 'trips' to home loc ####
### TOTAL ###
Choma5_clust_trips_shorter2_home <- Choma5_clust_trips_shorter2[which(Choma5_clust_trips_shorter2$home_clust2 == 1),]
Choma5_clust_trips_shorter2_home2 <- Choma5_clust_trips_shorter2_home[order(Choma5_clust_trips_shorter2_home$partid,
                                                                            Choma5_clust_trips_shorter2_home$clust,
                                                                            -Choma5_clust_trips_shorter2_home$trip_count_new),]
Choma5_clust_trips_shorter2_home3 <- Choma5_clust_trips_shorter2_home2[!duplicated(Choma5_clust_trips_shorter2_home2[,c(1,15)]),]

favstats(Choma5_clust_trips_shorter2_home3$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#   1  4      7 18.5  56 14.22034 15.03452 59       0

ggplot(data=Choma5_clust_trips_shorter2_home3, aes(x=trip_count_new)) +
  geom_histogram(binwidth = 5, fill = "red") +
  coord_cartesian(xlim=c(0,60))

### total time at home loc ####
Choma5_clust_trips_shorter2_home4 <- Choma5_clust_trips_shorter2_home[!duplicated(Choma5_clust_trips_shorter2_home$partid),]
favstats(Choma5_clust_trips_shorter2_home4$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
# 0.0156713 4.449881 20.42802 28.22408 116.9807 19.73786 17.44949 59       0

ggplot(data=Choma5_clust_trips_shorter2_home4, aes(x=(clust_trips_time_new))) +
  geom_histogram(binwidth = 5, fill = "red") +
  coord_cartesian(xlim=c(0,200))
#
### percent time at home ####
Choma5_clusts_shorter2_home$perc_time_home <- (Choma5_clusts_shorter2_home$clust_trips_time_new/Choma5_clusts_shorter2_home$part_trips_time_new)*100
favstats(Choma5_clusts_shorter2_home$perc_time_home)
#  min       Q1   median       Q3   max    mean       sd  n missing
#  0.08520987 13.64563 81.22594 95.40338 100 59.96667 39.03765 59       0
## without 3 people who spent 0% of time at home

Choma5_clusts_shorter2_home$partid[which(Choma5_clusts_shorter2_home$perc_time_home < 1)]
## still 3 people with <1% time at home
# "15084" "15085" "15097"

ggplot(data=Choma5_clusts_shorter2_home, aes(x=(perc_time_home))) +
  geom_histogram(binwidth = 10, fill = "black") +
  coord_cartesian(xlim=c(0,100))
##
# examine those with low percent time at home #####
parts_low_perc_home <- Choma5_clusts_shorter2_home$partid[which(Choma5_clusts_shorter2_home$perc_time_home < 15)]
temp_trips <- Choma5_clust_trips_shorter2[which(Choma5_clust_trips_shorter2$partid %in% parts_low_perc_home),]
temp_clusts <- Choma5_clusts_shorter2[which(Choma5_clusts_shorter2$partid %in% parts_low_perc_home),]
temp_clusts_no_home <- temp_clusts[which(temp_clusts$home_clust2 == 0),]

### number of locs
favstats(favstats(temp_clusts_no_home$clust~temp_clusts_no_home$partid)$n)
# min Q1  median    Q3 max     mean    sd    n  missing
#   7 12.25     18 20  22 16.14286 4.801557 14       0
# VS overall:
# min Q1 median Q3 max     mean       sd  n missing
#   1  7     12 15  17 10.77419 4.267373 62       0
hist(loc_counts2, xlim=c(0,25))
hist(favstats(temp_clusts_no_home$clust~temp_clusts_no_home$partid)$n, xlim=c(0,25), add=TRUE, col="red")
#
#### dist of locs
medians_km_all <- (favstats(Choma5_clusts_shorter2_no_home$hhdist_m_haversine~Choma5_clusts_shorter2_no_home$partid)$median)/1000
medians_km_low <- (favstats(temp_clusts_no_home$hhdist_m_haversine~temp_clusts_no_home$partid)$median)/1000
favstats(medians_km_low)
#    min    Q1 median    Q3   max  mean    sd  n missing
# 0.2448758 12.68832 14.69716 21.42752 22.63018 15.63178 6.084979 14       0

x <- seq(0,500, length.out=5000)
df <- data.frame(x = x, y1 = dnorm(x, mean(medians_km_all), sd(medians_km_all)),
                 y2 = dnorm(x, mean(medians_km_low), sd(medians_km_low)))
dists_km2 <- as.data.frame(cbind(c(medians_km_all, medians_km_low),
                                 c(rep("all",length(medians_km_all)),
                                   rep("low perc",length(medians_km_low)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) + 
  scale_fill_manual(values = c("grey","red")) +
  scale_color_manual(values = c("grey","red")) +
  coord_cartesian(xlim=c(0,30)) +
  scale_x_continuous(breaks=c(seq(0,30,5)))
##
### number of trips to locs ##
temp_trips_no_home3 <- Choma5_clust_trips_shorter2_no_home3[which(Choma5_clust_trips_shorter2_no_home3$partid %in% parts_low_perc_home),]
trips_part_median_all <- favstats(Choma5_clust_trips_shorter2_no_home3$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3$partid)$median
trips_part_median_low <- favstats(temp_trips_no_home3$trip_count_new ~ temp_trips_no_home3$partid)$median
favstats(trips_part_median_low)
# min  Q1 median Q3 max  mean    sd  n missing
#   1 1.625      2  2   2 1.75 0.4274252 14       0

x <- seq(1,500, length.out=5000)
df <- data.frame(x = x, y1 = dnorm(x, mean(trips_part_median_all), sd(trips_part_median_all)),
                 y2 = dnorm(x, mean(trips_part_median_low), sd(trips_part_median_low)))
dists_km2 <- as.data.frame(cbind(c(trips_part_median_all, trips_part_median_low),
                                 c(rep("all",length(trips_part_median_all)),
                                   rep("low_perc",length(trips_part_median_low)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("grey","red")) +
  scale_color_manual(values = c("grey","red")) +
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks=c(seq(0,10,1)))
##
### time per trip ##
temp_trips_no_home <- Choma5_clust_trips_shorter2_no_home[which(Choma5_clust_trips_shorter2_no_home$partid %in% parts_low_perc_home),]
trip_time_part_median_all <- favstats(Choma5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home$partid)$median
trip_time_part_median_low <- favstats(temp_trips_no_home$trip_time_new_hour ~ temp_trips_no_home$partid)$median
favstats(trip_time_part_median_low)
# min    Q1 median    Q3   max  mean    sd  n missing
# 1.040833 1.261042 1.589583 2.186493 3.834722 1.815059 0.7417849 14       0

x <- seq(0,500, length.out=5000)
df <- data.frame(x = x, y1 = dnorm(x, mean(trip_time_part_median_all), sd(trip_time_part_median_all)),
                 y2 = dnorm(x, mean(trip_time_part_median_low), sd(trip_time_part_median_low)))
dists_km2 <- as.data.frame(cbind(c(trip_time_part_median_all, trip_time_part_median_low),
                                 c(rep("all",length(trip_time_part_median_all)),
                                   rep("low_perc",length(trip_time_part_median_low)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) +
  scale_fill_manual(values = c("grey","red")) +
  scale_color_manual(values = c("grey","red")) +
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks=c(seq(0,10,1)))
##
### Average Direction of Travel #####
Choma5_clusts_shorter2_no_home$dir_from_home_bearing <- numeric(length=length(Choma5_clusts_shorter2_no_home$partid))

parts <- levels(as.factor(Choma5_clusts_shorter2_no_home$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Choma5_clusts_shorter2_no_home[which(Choma5_clusts_shorter2_no_home$partid == part),]
  for(j in 1:nrow(temp)){
    temp$dir_from_home_deg[j] <- bearing(c(temp$hh_long[j], temp$hh_lat[j]), c(temp$log_long[j], temp$log_lat[j]))
    if(temp$dir_from_home_deg[j] < 0){
      temp$dir_from_home_deg[j] <- 360+temp$dir_from_home_deg[j]
    }
  }
  Choma5_clusts_shorter2_no_home$dir_from_home_bearing[which(Choma5_clusts_shorter2_no_home$partid == part)] <- temp$dir_from_home_deg
}
## to get point on unit circle, must have angle theta from
## positive x axis, going positive to the counterclockwise direction
## these are in degrees from positive y going clockwise (bearings)
Choma5_clusts_shorter2_no_home$dir_from_home_deg <- 90 - Choma5_clusts_shorter2_no_home$dir_from_home_bearing
Choma5_clusts_shorter2_no_home$dir_from_home_deg[which(Choma5_clusts_shorter2_no_home$dir_from_home_deg < 0)] <- Choma5_clusts_shorter2_no_home$dir_from_home_deg[which(Choma5_clusts_shorter2_no_home$dir_from_home_deg < 0)] + 360
Choma5_clusts_shorter2_no_home$dir_from_home_rad <- deg2rad(Choma5_clusts_shorter2_no_home$dir_from_home_deg)

## average degree and magnitude
Choma5_clusts_shorter2_no_home$avg_dir_from_home_deg <- as.numeric(NA)
Choma5_clusts_shorter2_no_home$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_no_home$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Choma5_clusts_shorter2_no_home[which(Choma5_clusts_shorter2_no_home$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Choma5_clusts_shorter2_no_home$avg_dir_from_home_deg[which(Choma5_clusts_shorter2_no_home$partid == part)] <- avg_deg
  Choma5_clusts_shorter2_no_home$avg_dir_from_home_mag[which(Choma5_clusts_shorter2_no_home$partid == part)] <- avg_vec_dist/1000
}

Choma5_clusts_shorter2_no_home$avg_dir_from_home_bearing <- 90 - Choma5_clusts_shorter2_no_home$avg_dir_from_home_deg
Choma5_clusts_shorter2_no_home$avg_dir_from_home_bearing[which(Choma5_clusts_shorter2_no_home$avg_dir_from_home_bearing < 0)] <- Choma5_clusts_shorter2_no_home$avg_dir_from_home_bearing[which(Choma5_clusts_shorter2_no_home$avg_dir_from_home_bearing < 0)] + 360


temp <- Choma5_clusts_shorter2_no_home[!duplicated(Choma5_clusts_shorter2_no_home$partid),]

summary(temp$avg_dir_from_home_bearing)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.438  67.166 115.413 145.344 238.923 344.598 
summary(temp$avg_dir_from_home_mag)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1915  1.0580  6.4661  8.7563 14.4481 33.0336 

ggplot(data=temp, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(fill="black", col="black", position=position_dodge2(preserve = "single"), width=1) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                             "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw()
#
# ggplot(data=temp, aes(x=avg_dir_from_home_bearing)) +
#   geom_histogram(binwidth = 12.25) +
#   coord_polar(theta="x",start=0) +
#   scale_x_continuous(limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
#                                                                              "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW"))
#

## average degree only, not accounting for distance from home
Choma5_clusts_shorter2_no_home$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_no_home$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Choma5_clusts_shorter2_no_home[which(Choma5_clusts_shorter2_no_home$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Choma5_clusts_shorter2_no_home$avg_dir_from_home_deg2[which(Choma5_clusts_shorter2_no_home$partid == part)] <- avg_deg
}

Choma5_clusts_shorter2_no_home$avg_dir_from_home_values <- 90 - Choma5_clusts_shorter2_no_home$avg_dir_from_home_deg2
Choma5_clusts_shorter2_no_home$avg_dir_from_home_values[which(Choma5_clusts_shorter2_no_home$avg_dir_from_home_values < 0)] <- Choma5_clusts_shorter2_no_home$avg_dir_from_home_values[which(Choma5_clusts_shorter2_no_home$avg_dir_from_home_values < 0)] + 360
temp <- Choma5_clusts_shorter2_no_home[!duplicated(Choma5_clusts_shorter2_no_home$partid),]

median(temp$avg_dir_from_home_values)

temp$avg_dir_from_home_values[which(temp$avg_dir_from_home_values > 348.75)] <- -(360-temp$avg_dir_from_home_values[which(temp$avg_dir_from_home_values > 348.75)])
ggplot(data=temp, aes(x=avg_dir_from_home_values, y=..count..)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw()

##
### time spent at clusters vs elsewhere (travel, less important places, etc) ####
Choma5_clusts_shorter2_temp <- Choma5_clusts_shorter2[!duplicated(Choma5_clusts_shorter2$partid),]
Choma5_clusts_shorter2_temp <- Choma5_clusts_shorter2_temp[which(Choma5_clusts_shorter2_temp$partid != "30240"),]
Choma5_clusts_shorter2_temp$perc_time_in_all_clusts <- Choma5_clusts_shorter2_temp$part_trips_time_new/Choma5_clusts_shorter2_temp$total_time

# percent time in clusters
summary(Choma5_clusts_shorter2_temp$perc_time_in_all_clusts)
#    Min. 1st Qu.  Median    Mean  3rd Qu.    Max.
# 0.3665  0.8424  0.9335  0.8939  0.9807  1.0002 

# time outside clusters
summary(Choma5_clusts_shorter2_temp$total_time - Choma5_clusts_shorter2_temp$part_trips_time_new)
#       Min.  1st Qu.  Median      Mean    3rd Qu.  Max. 
# -0.00625  0.79480  2.09864  3.35031  4.63720 20.71190 
#
### Time at locations ####
favstats(Choma5_clusts_shorter2_no_home$clust_trips_time_new*24)
# min       Q1   median       Q3      max     mean       sd   n missing
# 0.2461111 1.084861 1.92125 4.797639 885.6553 27.36607 114.5253 738       0

clusts_median_Choma <- favstats((Choma5_clusts_shorter2_no_home$clust_trips_time_new*24)~Choma5_clusts_shorter2_no_home$partid)$median
favstats(clusts_median_Choma)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.7801389 1.369722 1.876111 2.4 7.662778 2.152006 1.260818 61       0

x <- seq(0,1000, length.out=1000)
df <- with(Choma5_clusts_shorter2_no_home, data.frame(x = x, y1 = dnorm(x, mean(clust_trips_time_new*24), sd(clust_trips_time_new*24)),
                                                           y2 = dnorm(x, mean(clusts_median_Choma), sd(clusts_median_Choma))))

clust_times_part2 <- as.data.frame(cbind(c(Choma5_clusts_shorter2_no_home$clust_trips_time_new*24, clusts_median_Choma),
                                         c(rep("overall",nrow(Choma5_clusts_shorter2_no_home)),
                                           rep("part_median",length(clusts_median_Choma)))))
colnames(clust_times_part2) <- c("values","type")
clust_times_part2$values <- as.numeric(as.character(clust_times_part2$values))
clust_times_part2$type <- factor(clust_times_part2$type)

ggplot(data=clust_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  #geom_line(data=df, aes(x=x, y=y1), color = "black") +
  #geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,900))


## just medians
ggplot(data=clust_times_part2[which(clust_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 1, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks=c(seq(0,10,1)))
#
#
####### Biting time metrics #####
##### at home ####
length(unique(Choma5_clusts_shorter2$partid)) # 62
length(unique(Choma5_clusts_shorter2$partid[which(!(is.na(Choma5_clusts_shorter2$clust_trips_time_new_biting_time)))])) # 62
## all parts spend time somewhere during biting hours
nrow(Choma5_clusts_shorter2_home) # 59
Choma5_clusts_shorter2_home_biting_places <- Choma5_clusts_shorter2_home[which(!(is.na(Choma5_clusts_shorter2_home$clust_trips_time_new_biting_time))),]
nrow(Choma5_clusts_shorter2_home_biting_places) # 53
## of 59 with home location, 6 don't spend any time there during biting hours (or unrecorded)
  ## for those 6, look at time spent at home in general

favstats(Choma5_clusts_shorter2_home_biting_places$clust_trips_time_new)
# min          Q1 median          Q3         max        mean          sd  n missing
# 0.7747916667 14.40002315  23.56 28.90712963 116.0581597 22.12279132 17.39327337 53       0
favstats(Choma5_clusts_shorter2_home_biting_places$clust_trips_time_new_biting_time)
# min          Q1      median          Q3         max        mean          sd  n missing
# 0.375 1.196689815 7.326423611 9.169236111 12.25751157 5.762152123 4.025270399 53       0

# percent of time at home that is during biting hours
favstats((Choma5_clusts_shorter2_home_biting_places$clust_trips_time_new_biting_time/Choma5_clusts_shorter2_home_biting_places$clust_trips_time_new)*100)
# min          Q1      median          Q3         max        mean          sd  n missing
# 3.052822634 16.49759339 28.09191842 35.92854984 51.72686804 27.29437263 12.51526021 53       0

##
##### outside home ####
### number locations ####
Choma5_clusts_shorter2_no_home_biting_places <- Choma5_clusts_shorter2_no_home[which(!(is.na(Choma5_clusts_shorter2_no_home$clust_trips_time_new_biting_time))),]
## only when >30 minutes during biting time
Choma5_clusts_shorter2_no_home_biting_places2 <- Choma5_clusts_shorter2_no_home_biting_places[which((Choma5_clusts_shorter2_no_home_biting_places$clust_trips_time_new_biting_time*24) >= 0.5),]
  ## removes 30 places, 3 participants

loc_counts <- favstats(Choma5_clusts_shorter2_no_home_biting_places2$clust~Choma5_clusts_shorter2_no_home_biting_places2$partid)$n
only_home <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2$partid))
loc_counts2 <- c(loc_counts,rep(0, only_home))
favstats(loc_counts2)
# min Q1 median Q3 max     mean       sd  n missing
#    0  0    1.5  3   6 1.919 1.777 62       0
summary(as.factor(loc_counts2))
#   0  1  2  3  4  5  6  
#  17 14  9 11  4  4  3 
## all those with 0 have biting time at home location (not just nowhere)
hist(loc_counts2, xlim=c(-1,10), breaks=c(seq(-1,10,by=1)))
##
### time per location ####
favstats((Choma5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5142 1.492  5.304 30.88 351.6 49.15 83.54 119       0

clusts_median_Choma <- (favstats(Choma5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2$partid)$median)*24
favstats(clusts_median_Choma)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.5928 2.621  8.546 91.95 291.3 53.11 76.77 45       0

x <- seq(0,1000, length.out=1000)
df <- with(Choma5_clusts_shorter2_no_home_biting_places2, data.frame(x = x, y1 = dnorm(x, mean(clust_trips_time_new_biting_time*24), sd(clust_trips_time_new_biting_time*24)),
                                                      y2 = dnorm(x, mean(clusts_median_Choma), sd(clusts_median_Choma))))

clust_times_part2 <- as.data.frame(cbind(c(Choma5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time*24, clusts_median_Choma),
                                         c(rep("overall",nrow(Choma5_clusts_shorter2_no_home_biting_places2)),
                                           rep("part_median",length(clusts_median_Choma)))))
colnames(clust_times_part2) <- c("values","type")
clust_times_part2$values <- as.numeric(as.character(clust_times_part2$values))
clust_times_part2$type <- factor(clust_times_part2$type)

ggplot(data=clust_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,360))

## just medians
ggplot(data=clust_times_part2[which(clust_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 10, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,300)) +
  scale_x_continuous(breaks=c(seq(0,300,20)))
#
#
### number trips ####
Choma5_clust_trips_shorter2_no_home_biting_places <- Choma5_clust_trips_shorter2_no_home[which(!(is.na(Choma5_clust_trips_shorter2_no_home$trip_time_new2_biting_time))),]
## only when >30 minutes during biting time
Choma5_clust_trips_shorter2_no_home_biting_places2 <- Choma5_clust_trips_shorter2_no_home_biting_places[which((Choma5_clust_trips_shorter2_no_home_biting_places$trip_time_new2_biting_time*24) >= 0.5),]
## removes 63 trips, 3 participants


Choma5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2 <- numeric(length=nrow(Choma5_clust_trips_shorter2_no_home_biting_places2))
parts <- levels(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places2$clust[which(Choma5_clust_trips_shorter2_no_home_biting_places2$partid == part)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    counts <- Choma5_clust_trips_shorter2_no_home_biting_places2$trip_count_new[which(Choma5_clust_trips_shorter2_no_home_biting_places2$partid == part & Choma5_clust_trips_shorter2_no_home_biting_places2$clust == clust)]
    counts <- as.factor(counts)
    levels(counts) <- c(1:nlevels(counts))
    Choma5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2[which(Choma5_clust_trips_shorter2_no_home_biting_places2$partid == part & Choma5_clust_trips_shorter2_no_home_biting_places2$clust == clust)] <- as.numeric(as.character(counts))
  }
}

Choma5_clust_trips_shorter2_no_home_biting_places2b <- Choma5_clust_trips_shorter2_no_home_biting_places2[order(Choma5_clust_trips_shorter2_no_home_biting_places2$partid,
                                                                                                                Choma5_clust_trips_shorter2_no_home_biting_places2$clust,
                                                                                                                -Choma5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2),]
Choma5_clust_trips_shorter2_no_home_biting_places3 <- Choma5_clust_trips_shorter2_no_home_biting_places2b[!duplicated(Choma5_clust_trips_shorter2_no_home_biting_places2b[,c(1,15)]),]

favstats(Choma5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  3  27 4.068 5.975 118       0
summary(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2))
#  1  2  3  4  5  6  7  8 10 12 13 14 15 17 21 22 23 27 
# 62 20  8  3  3  1  1  3  3  2  2  1  1  1  1  3  1  2 

### PER PART ###
trip_count_median <- favstats(Choma5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3$partid)$median
favstats(trip_count_median)
# min Q1 median Q3 max  mean  sd  n missing
#  1  1      2 3.5  27  3.6 5.147 45       0

### PLOT BOTH ###
x <- seq(0,500, length.out=1000)
df <- with(Choma5_clust_trips_shorter2_no_home_biting_places3, data.frame(x = x, y1 = dnorm(x, mean(trip_count_new2), sd(trip_count_new2)),
                                                            y2 = dnorm(x, mean(trip_count_median), sd(trip_count_median))))

trips_part2 <- as.data.frame(cbind(c(Choma5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2, trip_count_median),
                                   c(rep("overall",nrow(Choma5_clust_trips_shorter2_no_home_biting_places3)),
                                     rep("part_median",length(trip_count_median)))))
colnames(trips_part2) <- c("values","type")
trips_part2$values <- as.numeric(as.character(trips_part2$values))
trips_part2$type <- factor(trips_part2$type)

ggplot(data=trips_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,30))

## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,30))+ 
  scale_x_continuous(breaks=c(seq(0,30,2)))
#


### time per trip ####
favstats((Choma5_clust_trips_shorter2_no_home_biting_places2$trip_time_new2_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5006 2.921  8.031 10.02 141.5 12.17 17.5 480       0

clusts_median_Choma <- (favstats(Choma5_clust_trips_shorter2_no_home_biting_places2$trip_time_new2_biting_time~Choma5_clust_trips_shorter2_no_home_biting_places2$partid)$median)*24
favstats(clusts_median_Choma)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.8267 3.72  7.694 9.678 81.45 10.77 14.45 45       0

x <- seq(0,1000, length.out=1000)
df <- with(Choma5_clust_trips_shorter2_no_home_biting_places2, data.frame(x = x, y1 = dnorm(x, mean(trip_time_new2_biting_time*24), sd(trip_time_new2_biting_time*24)),
                                                                          y2 = dnorm(x, mean(clusts_median_Choma), sd(clusts_median_Choma))))

clust_times_part2 <- as.data.frame(cbind(c(Choma5_clust_trips_shorter2_no_home_biting_places2$trip_time_new2_biting_time*24, clusts_median_Choma),
                                         c(rep("overall",nrow(Choma5_clust_trips_shorter2_no_home_biting_places2)),
                                           rep("part_median",length(clusts_median_Choma)))))
colnames(clust_times_part2) <- c("values","type")
clust_times_part2$values <- as.numeric(as.character(clust_times_part2$values))
clust_times_part2$type <- factor(clust_times_part2$type)

ggplot(data=clust_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,140)) +
  scale_x_continuous(breaks=c(seq(0,200,20)))

## just medians
ggplot(data=clust_times_part2[which(clust_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 5, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,90)) +
  scale_x_continuous(breaks=c(seq(0,300,10)))
#
#

## just medians without 2 outliers
temp <- clust_times_part2[which(clust_times_part2$type == "part_median" & clust_times_part2$values < 40),]
ggplot(data=temp, aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 2, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,30)) +
  scale_x_continuous(breaks=c(seq(0,300,2)))
##### big group with ~7 hours --> 9pm-6am is 9 hours, so probably an overnight trip  



#
#
##### percent biting time spent at home vs elsewhere #####
Choma5_clusts_shorter2_biting_places <- Choma5_clusts_shorter2[which(!(is.na(Choma5_clusts_shorter2$clust_trips_time_new_biting_time))),]

Choma5_clusts_shorter2_biting_places$percent_home_biting_time <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_biting_places$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  percent <- Choma5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Choma5_clusts_shorter2_biting_places$partid == part & Choma5_clusts_shorter2_biting_places$home_clust2 == 1)]/
    sum(Choma5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Choma5_clusts_shorter2_biting_places$partid == part)])
  if(length(percent) == 0){
    percent <- 0
  }
  Choma5_clusts_shorter2_biting_places$percent_home_biting_time[which(Choma5_clusts_shorter2_biting_places$partid == part)] <- percent*100
}

Choma5_clusts_shorter2_biting_places2 <- Choma5_clusts_shorter2_biting_places[!duplicated(Choma5_clusts_shorter2_biting_places$partid),]
## 9 with 0%, but 3 of these don't have home location
homes <- Choma5_clusts_shorter2_home$partid
Choma5_clusts_shorter2_biting_places3 <- Choma5_clusts_shorter2_biting_places[which(Choma5_clusts_shorter2_biting_places$partid %in% homes),]
Choma5_clusts_shorter2_biting_places3b <- Choma5_clusts_shorter2_biting_places3[!duplicated(Choma5_clusts_shorter2_biting_places3$partid),]

favstats(Choma5_clusts_shorter2_biting_places3b$percent_home_biting_time)
favstats(Choma5_clusts_shorter2_biting_places3b$percent_home_biting_time[which(Choma5_clusts_shorter2_biting_places3b$percent_home_biting_time != 0)])

ggplot(data=Choma5_clusts_shorter2_biting_places3b, aes(x=(percent_home_biting_time))) +
  geom_histogram(binwidth = 10, fill = "black") +
  coord_cartesian(xlim=c(0,100)) +
  scale_x_continuous(breaks=c(seq(0,100,10)))
##
###############
###############
##### MUTASA #####
setwd("/Users/kschaber/Desktop/ICEMR_GPSLoggers/Data")
zimbabwe <- st_read("Zimbabwe.shp")
zimbabwe2 <- as_Spatial(zimbabwe)

point <- data.frame(lon=Mutasa5_clusts_shorter2$log_long, lat=Mutasa5_clusts_shorter2$log_lat)
sp2   <- SpatialPoints(point,proj4string=CRS(proj4string(zimbabwe2)))
Mutasa5_clusts_shorter2$in_zimbabwe <- as.factor(as.integer(gContains(zimbabwe2,sp2, byid=TRUE)))

moz_bound <- st_read("Mozambique.shp")
moz_bound2 <- as_Spatial(moz_bound)

point <- data.frame(lon=Mutasa5_clusts_shorter2$log_long, lat=Mutasa5_clusts_shorter2$log_lat)
sp2   <- SpatialPoints(point,proj4string=CRS(proj4string(moz_bound2)))
Mutasa5_clusts_shorter2$in_moz <- as.factor(as.integer(gContains(moz_bound2,sp2, byid=TRUE)))
setwd("/Users/kschaber/Desktop/ICEMR_GPSLoggers/Code/HDBSCAN/Outputs/New/3.17")


summary(Mutasa5_clusts_shorter2$in_zimbabwe) ## 55 clusters out of coutnry 
summary(Mutasa5_clusts_shorter2$in_moz) ## 52 clusters in Mozambique

Mutasa5_clusts_shorter2[which(Mutasa5_clusts_shorter2$in_zimbabwe == 0 & Mutasa5_clusts_shorter2$in_moz == 0),]
  ## 3 cluster locations just across the southern border in South Africa
      ## 1 person in each of their two measurements

Mutasa5_clusts_shorter2_moz <- Mutasa5_clusts_shorter2[which(Mutasa5_clusts_shorter2$in_moz == 1),]
sort(summary(as.factor(Mutasa5_clusts_shorter2_moz$partid)))
# 30182b 30184b 31064b  31068  33202  33311  33428  33429  33475  33570  33584  33591  33602  33603  33634  33683  33722 
#      3      3      4      4      4      1      6      8      3      2      5      1      1      1      1     13      3 

## 17 parts have clusters in Mozambique --> 5 have 1, 1 has 2, 4 have 3, 1 has 4, 2 have 5 and 2 have 8
favstats( Mutasa5_clusts_shorter2_moz$time_inclust ~Mutasa5_clusts_shorter2_moz$partid)
unique(Mutasa5_clusts_shorter2_moz$time_inclust)

#### OVERNIGHT SPECIFIC GPS (for Ellen) -- 10 parts ###
## 30182b --> once overnight (7/14-7/15)
## 30184b --> once overnight (7/14-7/15)
## 31064b --> over two nights (5/11-5/13 in Moz)
## 33202 --> over four nights (8/6-8/10 in Moz)
## 33311 --> over 2 nights (10/19-10/21 in Moz)
## 33428 --> over twenty-five nights (11/20-12/15 in Moz)
## 33429 --> over eighteen nights (11/27-12/15 in Moz)
## 33475 --> over thirteen nights (1/30-2/12 in Moz)
## 33584 --> over four nights (3/18-3/22 in Moz)
## 33683 --> over thirty-four nights (5/17-6/20 (goes offline and back on in Moz on 6/23, then leaves) in Moz)

## 7 with clusters in Moz, but don't visit overnight
unique(Mutasa5_clusts_shorter2_moz$partid)
# 31068 33570 33591 33602 33603 33634 33722
Mutasa5_clusts_shorter2_moz2 <- Mutasa5_clusts_shorter2_moz[which(Mutasa5_clusts_shorter2_moz$partid %in% c("31068", "33570", "33591", "33602", "33603", "33634", "33722")),]
summary(Mutasa5_clusts_shorter2_moz2$time_inclust)
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.06905 0.17032 0.18767 0.24844 0.28531 0.78356 
summary(Mutasa5_clusts_shorter2_moz2$hhdist_m_haversine)/1000
#   Min.  1st Qu.   Median     Mean   3rd Qu.    Max. 
# 1.575    3.109     3.935   13.281    6.629   86.905


##### in all of gps points, ## 39 parts of 228 parts total had GPS in Mozambique
## ## 27 with >5 points in Mozambique

##
#### Get values for biting time period ####
Mutasa5_new2 <- Mutasa5_new
Mutasa5_new2$biting_time <- factor(0, levels = c(0,1))
Mutasa5_new2$biting_time[which(Mutasa5_new2$time_new_hour %in% c(21:23,0:5))] <- 1

Mutasa5_new2$trip_time_new2_biting_time <- as.numeric(NA)
Mutasa5_new2$clust_trips_time_new_biting_time <- as.numeric(NA)

parts <- levels(as.factor(Mutasa5_new2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Mutasa5_new2$clust[which(Mutasa5_new2$partid == part)]))
  for(jj in 1:length(clusts)){
    clust <- clusts[jj]
    trips <- levels(as.factor(Mutasa5_new2$trip_count_new[which(Mutasa5_new2$partid == part & Mutasa5_new2$clust == clust)]))
    temp3 <- Mutasa5_new2[which(Mutasa5_new2$partid == part & Mutasa5_new2$clust == clust & Mutasa5_new2$biting_time == 1),]
    if(nrow(temp3) > 0){
      xx2 <- aggregate(time_atpoint_new ~ trip_count_new, temp3, sum)
      for(kk in 1:nrow(xx2)){
        Mutasa5_new2$trip_time_new2_biting_time[which(Mutasa5_new2$partid == part & 
                                                        Mutasa5_new2$clust == clust &
                                                        Mutasa5_new2$trip_count_new == xx2$trip_count_new[kk])] <- xx2$time_atpoint_new[kk]
      }
    }
  }
  temp4 <- Mutasa5_new2[which(Mutasa5_new2$partid == part),]
  temp4b <- temp4[!duplicated(temp4[c(1,22,26)]),]
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_biting_time))),])>0){
    biting <- aggregate(trip_time_new2_biting_time ~ clust, temp4b, sum)
    if(nrow(biting) > 0){
      for(kk in 1:nrow(biting)){
        Mutasa5_new2$clust_trips_time_new_biting_time[which(Mutasa5_new2$partid == part & Mutasa5_new2$clust == biting$clust[kk] & !(is.na(Mutasa5_new2$trip_count_new)))] <- biting$trip_time_new2_biting_time[kk]
      }
    }
  }
}
Mutasa5_new2$part_trips_time_new_biting_time <- as.numeric(0)
temp <- Mutasa5_new2[!duplicated(Mutasa5_new2[c(1,22)]),]
tot <- aggregate(clust_trips_time_new_biting_time ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Mutasa5_new2$part_trips_time_new_biting_time[which(Mutasa5_new2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_biting_time[ii]
}

prop.table(table(Mutasa5_new2$home_clust2, Mutasa5_new2$biting_time, deparse.level = 2), margin=2)*100
#                         Mutasa5_new2$biting_time
# Mutasa5_new2$home_clust2        0        1
# 0                        40.21593 19.01075
# 1                        59.78407 80.98925
prop.table(table(Mutasa5_new2$home_clust2, Mutasa5_new2$biting_time, deparse.level = 2), margin=1)*100
#                         Mutasa5_new2$biting_time
# Mutasa5_new2$home_clust2        0        1
# 0                        87.05616 12.94384
# 1                        70.12165 29.87835

#### 81% of observations taken during biting times were in the home cluster
#### 30% of obervations in the home cluster were taken during biting times

Mutasa5_clusts_shorter_temp <- Mutasa5_new2[order(Mutasa5_new2$partid, Mutasa5_new2$clust),]
Mutasa5_clusts_shorter_temp2 <- Mutasa5_clusts_shorter_temp[!duplicated(Mutasa5_clusts_shorter_temp[,c(1,22)]),c(1:7,16,17,12,13,19:25,29,30:35,46:50,41:45,51,55,54)]

Mutasa5_clust_trips_shorter_temp <- Mutasa5_new2[order(Mutasa5_new2$partid, Mutasa5_new2$clust, Mutasa5_new2$trip_count),]
Mutasa5_clust_trips_shorter_temp2 <- Mutasa5_clust_trips_shorter_temp[!duplicated(Mutasa5_clust_trips_shorter_temp[,c(1,22,26)]),c(1:7,16,17,12,13,19:22,24,26:28,29,30:35,46:50,41:45,36:40,51,55,54,53)]


Mutasa5_clusts_shorter2$part_trips_time_new_biting_time <- Mutasa5_clusts_shorter_temp2$part_trips_time_new_biting_time
Mutasa5_clusts_shorter2$clust_trips_time_new_biting_time <- Mutasa5_clusts_shorter_temp2$clust_trips_time_new_biting_time

Mutasa5_clust_trips_shorter2$part_trips_time_new_biting_time <- Mutasa5_clust_trips_shorter_temp2$part_trips_time_new_biting_time
Mutasa5_clust_trips_shorter2$clust_trips_time_new_biting_time <- Mutasa5_clust_trips_shorter_temp2$clust_trips_time_new_biting_time
Mutasa5_clust_trips_shorter2$trip_time_new2_biting_time <- Mutasa5_clust_trips_shorter_temp2$trip_time_new2_biting_time
##
### number of locations (without home) ####
Mutasa5_clusts_shorter2_no_home <- Mutasa5_clusts_shorter2[which(Mutasa5_clusts_shorter2$home_clust2 == 0),]
loc_counts <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)$n
only_home <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home$partid))
loc_counts2 <- c(loc_counts,rep(0, only_home))
favstats(loc_counts2)
# min Q1 median Q3 max     mean       sd  n missing
#  0  3.25      7 11  25 7.934066 5.695403 182       0
hist(loc_counts2, xlim=c(0,27), breaks=c(0,5,10,15,20,25))
#
### distance of locations (without home) ####
### TOTAL ###
Mutasa5_clusts_shorter2_no_home$hhdist_km_haversine <- Mutasa5_clusts_shorter2_no_home$hhdist_m_haversine/1000
favstats(Mutasa5_clusts_shorter2_no_home$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#  2.745604e-05 0.5001869 1.325419 4.292403 481.4493 16.92257 54.57856 1444       0

### PER PART ###
medians_km <- (favstats(Mutasa5_clusts_shorter2_no_home$hhdist_m_haversine~Mutasa5_clusts_shorter2_no_home$partid)$median)/1000
favstats(medians_km) # km
#          min        Q1   median       Q3      max     mean       sd  n missing
# 0.08078032 0.624203 1.200587 2.59482 299.3132 10.47779 37.43025 176       0

medians_km_short <- medians_km[-which(medians_km > 100)] ## 6 outliers
favstats(medians_km_short)
# min        Q1    median       Q3      max     mean       sd  n missing
# 0.08078032 0.6155971 1.073367 2.262333 57.42812 4.071488 9.115249 170       0

x <- seq(0,500, length.out=500)
df <- with(Mutasa5_clusts_shorter2_no_home, data.frame(x = x, y1 = dnorm(x, mean(hhdist_km_haversine), sd(hhdist_km_haversine)),
                                                       y2 = dnorm(x, mean(medians_km_short), sd(medians_km_short))))

dists_km2 <- as.data.frame(cbind(c(Mutasa5_clusts_shorter2_no_home$hhdist_km_haversine, medians_km_short),
                                 c(rep("overall",nrow(Mutasa5_clusts_shorter2_no_home)),
                                   rep("part_median",length(medians_km_short)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,500))
#

## just medians
ggplot(data=dists_km2[which(dists_km2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 5, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=5) +
  coord_cartesian(xlim=c(0,60))
#
### number of trips (without home) ####
### TOTAL ###
Mutasa5_clust_trips_shorter2_no_home <- Mutasa5_clust_trips_shorter2[which(Mutasa5_clust_trips_shorter2$home_clust2 == 0),]
Mutasa5_clust_trips_shorter2_no_home2 <- Mutasa5_clust_trips_shorter2_no_home[order(Mutasa5_clust_trips_shorter2_no_home$partid,
                                                                                    Mutasa5_clust_trips_shorter2_no_home$clust,
                                                                                    -Mutasa5_clust_trips_shorter2_no_home$trip_count_new),]
Mutasa5_clust_trips_shorter2_no_home3 <- Mutasa5_clust_trips_shorter2_no_home2[!duplicated(Mutasa5_clust_trips_shorter2_no_home2[,c(1,15)]),]

favstats(Mutasa5_clust_trips_shorter2_no_home3$trip_count_new)
# min Q1 median  Q3  max     mean       sd     n missing
# 1  1      2  3  43 2.581025 3.609635 1444       0

### PER PART ###
trips_part_median <- favstats(Mutasa5_clust_trips_shorter2_no_home3$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3$partid)$median
favstats(trips_part_median)
#    min Q1 median Q3 max     mean       sd  n missing
#  1  1    1.5  2  13 1.892045 1.662956 176       0

### PLOT BOTH ###
x <- seq(0,100, length.out=500)
df <- with(Mutasa5_clust_trips_shorter2_no_home3, data.frame(x = x, y1 = dnorm(x, mean(trip_count_new), sd(trip_count_new)),
                                                             y2 = dnorm(x, mean(trips_part_median), sd(trips_part_median))))

trips_part2 <- as.data.frame(cbind(c(Mutasa5_clust_trips_shorter2_no_home3$trip_count_new, trips_part_median),
                                   c(rep("overall",nrow(Mutasa5_clust_trips_shorter2_no_home3)),
                                     rep("part_median",length(trips_part_median)))))
colnames(trips_part2) <- c("values","type")
trips_part2$values <- as.numeric(as.character(trips_part2$values))
trips_part2$type <- factor(trips_part2$type)

ggplot(data=trips_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,35))


## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=2) +
  coord_cartesian(xlim=c(0,14))+ 
  scale_x_continuous(breaks=c(seq(0,14,2)))
#
### time per trip ####
### TOTAL ###
Mutasa5_clust_trips_shorter2_no_home$trip_time_new_hour <- Mutasa5_clust_trips_shorter2_no_home$trip_time_new2*24
favstats(Mutasa5_clust_trips_shorter2_no_home$trip_time_new_hour)
#         min    Q1   median       Q3      max     mean       sd    n missing
# 0.4413889 1.384583 2.613333 5.324861 1218.114 6.898847 34.61174 3727       0

### PER PART ###
trip_time_part_median <- favstats(Mutasa5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home$partid)$median
favstats(trip_time_part_median)
#       min       Q1   median       Q3      max     mean       sd  n missing
# 0.9602778 2.062882 2.843333 3.930139 1218.114 10.45214 91.60288 176       0

### PLOT BOTH ###
## REMOVE OUTLIER ##
trip_time_part_median_short <- trip_time_part_median[which(trip_time_part_median < 30)]
favstats(trip_time_part_median_short)
# min       Q1   median       Q3      max     mean       sd   n missing
# 0.9602778 2.047535 2.8275 3.866944 12.22667 3.376434 2.052917 174       0

x <- seq(0,5000, length.out=50000)
df <- with(Mutasa5_clust_trips_shorter2_no_home, data.frame(x = x, y1 = dnorm(x, mean(trip_time_new_hour), sd(trip_time_new_hour)),
                                                            y2 = dnorm(x, mean(trip_time_part_median_short), sd(trip_time_part_median_short))))

trip_times_part2 <- as.data.frame(cbind(c(Mutasa5_clust_trips_shorter2_no_home$trip_time_new_hour, trip_time_part_median_short),
                                        c(rep("overall",nrow(Mutasa5_clust_trips_shorter2_no_home)),
                                          rep("part_median",length(trip_time_part_median_short)))))
colnames(trip_times_part2) <- c("values","type")
trip_times_part2$values <- as.numeric(as.character(trip_times_part2$values))
trip_times_part2$type <- factor(trip_times_part2$type)

ggplot(data=trip_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "darkgrey") +
  coord_cartesian(xlim=c(0,1250))
#

## just medians
ggplot(data=trip_times_part2[which(trip_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 1, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1.5) +
  coord_cartesian(xlim=c(0,15))
#
### Home Loc ####
Mutasa5_clusts_shorter2_home <- Mutasa5_clusts_shorter2[which(Mutasa5_clusts_shorter2$home_clust2 == 1),]
no_home <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_home$partid))# 0
## none with no homes 
### number of 'trips' to home loc ####
### TOTAL ###
Mutasa5_clust_trips_shorter2_home <- Mutasa5_clust_trips_shorter2[which(Mutasa5_clust_trips_shorter2$home_clust2 == 1),]
Mutasa5_clust_trips_shorter2_home2 <- Mutasa5_clust_trips_shorter2_home[order(Mutasa5_clust_trips_shorter2_home$partid,
                                                                              Mutasa5_clust_trips_shorter2_home$clust,
                                                                              -Mutasa5_clust_trips_shorter2_home$trip_count_new),]
Mutasa5_clust_trips_shorter2_home3 <- Mutasa5_clust_trips_shorter2_home2[!duplicated(Mutasa5_clust_trips_shorter2_home2[,c(1,15)]),]

favstats(Mutasa5_clust_trips_shorter2_home3$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#  1  6     11 25  49 15.40659 11.91357 182       0

ggplot(data=Mutasa5_clust_trips_shorter2_home3, aes(x=trip_count_new)) +
  geom_histogram(binwidth = 5, fill = "red") +
  coord_cartesian(xlim=c(0,100))

### total time at home loc ####
Mutasa5_clust_trips_shorter2_home4 <- Mutasa5_clust_trips_shorter2_home[!duplicated(Mutasa5_clust_trips_shorter2_home$partid),]
favstats(Mutasa5_clust_trips_shorter2_home4$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
#  5.420035 18.90386 24.81499 30.44999 132.13 26.22769 13.32921 182       0
favstats(Mutasa5_clust_trips_shorter2_home4$clust_trips_time_new[which(Mutasa5_clust_trips_shorter2_home4$partid != "30240")])
#        min       Q1   median       Q3      max     mean       sd  n missing
# 5.420035 18.88848 24.74998 30.09 65.96001 25.6426 10.77049 181       0

ggplot(data=Mutasa5_clust_trips_shorter2_home4[which(Mutasa5_clust_trips_shorter2_home4$partid != "30240"),], aes(x=clust_trips_time_new)) +
  geom_histogram(binwidth = 5, fill = "red") +
  coord_cartesian(xlim=c(0,70))
#
### percent time at home ####
Mutasa5_clusts_shorter2_home$perc_time_home <- (Mutasa5_clusts_shorter2_home$clust_trips_time_new/Mutasa5_clusts_shorter2_home$part_trips_time_new)*100
favstats(Mutasa5_clusts_shorter2_home$perc_time_home)
#  min       Q1   median       Q3   max    mean       sd  n missing
# 14.3927 72.97264 88.41037 96.42304 100 81.78818 18.65673 182       0

Mutasa5_clusts_shorter2_home$partid[which(Mutasa5_clusts_shorter2_home$perc_time_home < 5)]
## none

ggplot(data=Mutasa5_clusts_shorter2_home, aes(x=(perc_time_home))) +
  geom_histogram(binwidth = 10, fill = "black") +
  coord_cartesian(xlim=c(0,100))
##
### Average Direction of Travel #####
Mutasa5_clusts_shorter2_no_home$dir_from_home_bearing <- numeric(length=length(Mutasa5_clusts_shorter2_no_home$partid))

parts <- levels(as.factor(Mutasa5_clusts_shorter2_no_home$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Mutasa5_clusts_shorter2_no_home[which(Mutasa5_clusts_shorter2_no_home$partid == part),]
  for(j in 1:nrow(temp)){
    temp$dir_from_home_deg[j] <- bearing(c(temp$hh_long[j], temp$hh_lat[j]), c(temp$log_long[j], temp$log_lat[j]))
    if(temp$dir_from_home_deg[j] < 0){
      temp$dir_from_home_deg[j] <- 360+temp$dir_from_home_deg[j]
    }
  }
  Mutasa5_clusts_shorter2_no_home$dir_from_home_bearing[which(Mutasa5_clusts_shorter2_no_home$partid == part)] <- temp$dir_from_home_deg
}
## to get point on unit circle, must have angle theta from
## positive x axis, going positive to the counterclockwise direction
## these are in degrees from positive y going clockwise (bearings)
Mutasa5_clusts_shorter2_no_home$dir_from_home_deg <- 90 - Mutasa5_clusts_shorter2_no_home$dir_from_home_bearing
Mutasa5_clusts_shorter2_no_home$dir_from_home_deg[which(Mutasa5_clusts_shorter2_no_home$dir_from_home_deg < 0)] <- Mutasa5_clusts_shorter2_no_home$dir_from_home_deg[which(Mutasa5_clusts_shorter2_no_home$dir_from_home_deg < 0)] + 360
Mutasa5_clusts_shorter2_no_home$dir_from_home_rad <- deg2rad(Mutasa5_clusts_shorter2_no_home$dir_from_home_deg)

## average degree and magnitude
Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_deg <- as.numeric(NA)
Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_no_home$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Mutasa5_clusts_shorter2_no_home[which(Mutasa5_clusts_shorter2_no_home$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_deg[which(Mutasa5_clusts_shorter2_no_home$partid == part)] <- avg_deg
  Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_mag[which(Mutasa5_clusts_shorter2_no_home$partid == part)] <- avg_vec_dist/1000
}

Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_bearing <- 90 - Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_deg
Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_bearing[which(Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_bearing < 0)] <- Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_bearing[which(Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_bearing < 0)] + 360


temp <- Mutasa5_clusts_shorter2_no_home[!duplicated(Mutasa5_clusts_shorter2_no_home$partid),]

ggplot(data=temp, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(fill="black", col="black", position=position_dodge2(preserve = "single"), width=1) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw()
#

summary(temp$avg_dir_from_home_bearing)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.064 114.562 219.978 200.400 283.178 359.737 
summary(temp$avg_dir_from_home_mag)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.03371   0.60786   1.42987  13.38188   5.65708 210.00803 

temp2 <- temp[which(temp$avg_dir_from_home_mag < 75),] ## 164 of 176 parts

summary(temp2$avg_dir_from_home_bearing)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.064 108.505 213.498 196.622 275.108 359.737 
summary(temp2$avg_dir_from_home_mag)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03371  0.56605  1.26050  5.19688  3.52588 65.41050 

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(fill="black", col="black", position=position_dodge2(preserve = "single"), width=1) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw()
#



## average degree only, not accounting for distance from home
Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_no_home$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Mutasa5_clusts_shorter2_no_home[which(Mutasa5_clusts_shorter2_no_home$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_deg2[which(Mutasa5_clusts_shorter2_no_home$partid == part)] <- avg_deg
}

Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_values <- 90 - Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_deg2
Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_values[which(Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_values < 0)] <- Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_values[which(Mutasa5_clusts_shorter2_no_home$avg_dir_from_home_values < 0)] + 360
temp <- Mutasa5_clusts_shorter2_no_home[!duplicated(Mutasa5_clusts_shorter2_no_home$partid),]
median(temp$avg_dir_from_home_values)

temp$avg_dir_from_home_values[which(temp$avg_dir_from_home_values > 348.75)] <- -(360-temp$avg_dir_from_home_values[which(temp$avg_dir_from_home_values > 348.75)])
ggplot(data=temp, aes(x=avg_dir_from_home_values, y=..count..)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw()

####
### time spent at clusters vs elsewhere (travel, less important places, etc) ####
Mutasa5_clusts_shorter2_temp <- Mutasa5_clusts_shorter2[!duplicated(Mutasa5_clusts_shorter2$partid),]
Mutasa5_clusts_shorter2_temp <- Mutasa5_clusts_shorter2_temp[which(Mutasa5_clusts_shorter2_temp$partid != "30240"),]
Mutasa5_clusts_shorter2_temp$perc_time_in_all_clusts <- Mutasa5_clusts_shorter2_temp$part_trips_time_new/Mutasa5_clusts_shorter2_temp$total_time

# percent time in clusters
summary(Mutasa5_clusts_shorter2_temp$perc_time_in_all_clusts)
#    Min. 1st Qu.  Median    Mean  3rd Qu.    Max.
# 0.3665  0.8424  0.9335  0.8939  0.9807  1.0002 

# time outside clusters
summary(Mutasa5_clusts_shorter2_temp$total_time - Mutasa5_clusts_shorter2_temp$part_trips_time_new)
#       Min.  1st Qu.  Median      Mean    3rd Qu.  Max. 
# -0.00625  0.79480  2.09864  3.35031  4.63720 20.71190 

## part 30320b --> ~1/3 of points considered 'outliers' --> seem to be around home, but spread out enough to not count
#

### Time at locations ####
favstats(Mutasa5_clusts_shorter2_no_home$clust_trips_time_new*24)
# min       Q1   median       Q3      max    mean       sd    n missing
# 0.4802778 2.279583 4.154167 9.009444 1537.291 17.8061 73.87933 1444       0
clusts_median_Mutasa <- favstats((Mutasa5_clusts_shorter2_no_home$clust_trips_time_new*24)~Mutasa5_clusts_shorter2_no_home$partid)$median
favstats(clusts_median_Mutasa)
# min       Q1   median       Q3      max     mean      sd   n missing
# 1.433611 3.263889 4.202153 6.139792 1218.114 15.29516 93.3154 176       0
clusts_median_Mutasa2 <- clusts_median_Mutasa[which(clusts_median_Mutasa < 150)]
favstats(clusts_median_Mutasa2)
# min       Q1  median       Q3      max     mean       sd   n missing
# 1.433611 3.251667 4.19125 6.007118 110.4386 7.347391 13.94902 174       0

x <- seq(0,1000, length.out=1000)
df <- with(Mutasa5_clusts_shorter2_no_home, data.frame(x = x, y1 = dnorm(x, mean(clust_trips_time_new*24), sd(clust_trips_time_new*24)),
                                                      y2 = dnorm(x, mean(clusts_median_Mutasa2), sd(clusts_median_Mutasa2))))

clust_times_part2 <- as.data.frame(cbind(c(Mutasa5_clusts_shorter2_no_home$clust_trips_time_new*24, clusts_median_Mutasa2),
                                         c(rep("overall",nrow(Mutasa5_clusts_shorter2_no_home)),
                                           rep("part_median",length(clusts_median_Mutasa2)))))
colnames(clust_times_part2) <- c("values","type")
clust_times_part2$values <- as.numeric(as.character(clust_times_part2$values))
clust_times_part2$type <- factor(clust_times_part2$type)

ggplot(data=clust_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  #geom_line(data=df, aes(x=x, y=y1), color = "black") +
  #geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,1600))

## just medians
ggplot(data=clust_times_part2[which(clust_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 10, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=7) +
  coord_cartesian(xlim=c(0,120)) +
  scale_x_continuous(breaks=c(seq(0,200,10)))
#
#
####### Biting time metrics #####
##### at home ####
length(unique(Mutasa5_clusts_shorter2$partid)) # 182
length(unique(Mutasa5_clusts_shorter2$partid[which(!(is.na(Mutasa5_clusts_shorter2$clust_trips_time_new_biting_time)))])) # 182
## all parts spend time somewhere during biting hours
nrow(Mutasa5_clusts_shorter2_home) # 182
Mutasa5_clusts_shorter2_home_biting_places <- Mutasa5_clusts_shorter2_home[which(!(is.na(Mutasa5_clusts_shorter2_home$clust_trips_time_new_biting_time))),]
nrow(Mutasa5_clusts_shorter2_home_biting_places) # 53
## of 182 with home location, all spend time there during biting hours

favstats(Mutasa5_clusts_shorter2_home_biting_places$clust_trips_time_new)
# min          Q1 median          Q3         max        mean          sd  n missing
# 5.416 18.91  24.81 30.47 83.36 25.95 11.56 182       0
favstats(Mutasa5_clusts_shorter2_home_biting_places$clust_trips_time_new_biting_time)
# min          Q1      median          Q3         max        mean          sd  n missing
# 0.9424 5.783  7.966 9.772 23.81 8.053 3.558 182       0

# percent of time at home that is during biting hours
favstats((Mutasa5_clusts_shorter2_home_biting_places$clust_trips_time_new_biting_time/Mutasa5_clusts_shorter2_home_biting_places$clust_trips_time_new)*100)
# min          Q1      median          Q3         max        mean          sd  n missing
# 2.79 26.08  33.33 39.07 62.47 32.96 11.1 182       0
##
##### outside home ####
### number locations ####
Mutasa5_clusts_shorter2_no_home_biting_places <- Mutasa5_clusts_shorter2_no_home[which(!(is.na(Mutasa5_clusts_shorter2_no_home$clust_trips_time_new_biting_time))),]
## only when >30 minutes during biting time
Mutasa5_clusts_shorter2_no_home_biting_places2 <- Mutasa5_clusts_shorter2_no_home_biting_places[which((Mutasa5_clusts_shorter2_no_home_biting_places$clust_trips_time_new_biting_time*24) >= 0.5),]
## removes 30 places, 3 participants

loc_counts <- favstats(Mutasa5_clusts_shorter2_no_home_biting_places2$clust~Mutasa5_clusts_shorter2_no_home_biting_places2$partid)$n
only_home <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2$partid))
loc_counts2 <- c(loc_counts,rep(0, only_home))
favstats(loc_counts2)
# min Q1 median Q3 max     mean       sd  n missing
#    0  0      1  2   8 1.39 1.764 182       0
summary(as.factor(loc_counts2))
#   0  1  2  3  4  5  6  7  8 
#  74 43 32 14  9  1  4  1  4 
## all those with 0 have biting time at home location (not just nowhere)
hist(loc_counts2, xlim=c(-1,10), breaks=c(seq(-1,10,by=1)))


### time per location ####
favstats((Mutasa5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5072 1.735  7.418 15.57 504.1 21.25 47.85 253       0

clusts_median_Mutasa <- (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2$partid)$median)*24
favstats(clusts_median_Mutasa)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.5831 3.81  8.444 17.41 504.1 24.53 58.43 108       0

x <- seq(0,1000, length.out=1000)
df <- with(Mutasa5_clusts_shorter2_no_home_biting_places2, data.frame(x = x, y1 = dnorm(x, mean(clust_trips_time_new_biting_time*24), sd(clust_trips_time_new_biting_time*24)),
                                                                     y2 = dnorm(x, mean(clusts_median_Mutasa), sd(clusts_median_Mutasa))))

clust_times_part2 <- as.data.frame(cbind(c(Mutasa5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time*24, clusts_median_Mutasa),
                                         c(rep("overall",nrow(Mutasa5_clusts_shorter2_no_home_biting_places2)),
                                           rep("part_median",length(clusts_median_Mutasa)))))
colnames(clust_times_part2) <- c("values","type")
clust_times_part2$values <- as.numeric(as.character(clust_times_part2$values))
clust_times_part2$type <- factor(clust_times_part2$type)

ggplot(data=clust_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,525))

## just medians
ggplot(data=clust_times_part2[which(clust_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 10, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,520)) +
  scale_x_continuous(breaks=c(seq(0,520,100)))
#
#
## just medians without 1 outlier
temp <- clust_times_part2[which(clust_times_part2$type == "part_median" & clust_times_part2$values <400),]
ggplot(data=temp, aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 10, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1.5) +
  coord_cartesian(xlim=c(0,234)) +
  scale_x_continuous(breaks=c(seq(0,520,20)))
#
#

### number trips ####
Mutasa5_clust_trips_shorter2_no_home_biting_places <- Mutasa5_clust_trips_shorter2_no_home[which(!(is.na(Mutasa5_clust_trips_shorter2_no_home$trip_time_new2_biting_time))),]
## only when >30 minutes during biting time
Mutasa5_clust_trips_shorter2_no_home_biting_places2 <- Mutasa5_clust_trips_shorter2_no_home_biting_places[which((Mutasa5_clust_trips_shorter2_no_home_biting_places$trip_time_new2_biting_time*24) >= 0.5),]
## removes 63 trips, 3 participants


Mutasa5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2 <- numeric(length=nrow(Mutasa5_clust_trips_shorter2_no_home_biting_places2))
parts <- levels(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places2$clust[which(Mutasa5_clust_trips_shorter2_no_home_biting_places2$partid == part)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    counts <- Mutasa5_clust_trips_shorter2_no_home_biting_places2$trip_count_new[which(Mutasa5_clust_trips_shorter2_no_home_biting_places2$partid == part & Mutasa5_clust_trips_shorter2_no_home_biting_places2$clust == clust)]
    counts <- as.factor(counts)
    levels(counts) <- c(1:nlevels(counts))
    Mutasa5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2[which(Mutasa5_clust_trips_shorter2_no_home_biting_places2$partid == part & Mutasa5_clust_trips_shorter2_no_home_biting_places2$clust == clust)] <- as.numeric(as.character(counts))
  }
}

Mutasa5_clust_trips_shorter2_no_home_biting_places2b <- Mutasa5_clust_trips_shorter2_no_home_biting_places2[order(Mutasa5_clust_trips_shorter2_no_home_biting_places2$partid,
                                                                                                                Mutasa5_clust_trips_shorter2_no_home_biting_places2$clust,
                                                                                                                -Mutasa5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2),]
Mutasa5_clust_trips_shorter2_no_home_biting_places3 <- Mutasa5_clust_trips_shorter2_no_home_biting_places2b[!duplicated(Mutasa5_clust_trips_shorter2_no_home_biting_places2b[,c(1,15)]),]

favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#  1  1      1  2  22 2.19 3.14 253       0
summary(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2))
#   1   2   3   4   5   6   7   8   9  10  11  13  14  22 
# 175  36  13   6   1   3   2   4   4   1   2   1   2   3

### PER PART ###
trip_count_median <- favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3$partid)$median
favstats(trip_count_median)
# min Q1 median Q3 max  mean  sd  n missing
# 1  1      1  2  13 1.801 1.969 108       0

### PLOT BOTH ###
x <- seq(0,500, length.out=1000)
df <- with(Mutasa5_clust_trips_shorter2_no_home_biting_places3, data.frame(x = x, y1 = dnorm(x, mean(trip_count_new2), sd(trip_count_new2)),
                                                                          y2 = dnorm(x, mean(trip_count_median), sd(trip_count_median))))

trips_part2 <- as.data.frame(cbind(c(Mutasa5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2, trip_count_median),
                                   c(rep("overall",nrow(Mutasa5_clust_trips_shorter2_no_home_biting_places3)),
                                     rep("part_median",length(trip_count_median)))))
colnames(trips_part2) <- c("values","type")
trips_part2$values <- as.numeric(as.character(trips_part2$values))
trips_part2$type <- factor(trips_part2$type)

ggplot(data=trips_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,25))

## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1.5) +
  coord_cartesian(xlim=c(0,14))+ 
  scale_x_continuous(breaks=c(seq(0,15,2)))
#


### time per trip ####
favstats((Mutasa5_clust_trips_shorter2_no_home_biting_places2$trip_time_new2_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5058 1.546  6.218 8.829 266.1 9.687 21.86 554       0

clusts_median_Mutasa <- (favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places2$trip_time_new2_biting_time~Mutasa5_clust_trips_shorter2_no_home_biting_places2$partid)$median)*24
favstats(clusts_median_Mutasa)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.5831 2.24  7.445 8.986 234.9 13.18 29.81 108       0

x <- seq(0,1000, length.out=1000)
df <- with(Mutasa5_clust_trips_shorter2_no_home_biting_places2, data.frame(x = x, y1 = dnorm(x, mean(trip_time_new2_biting_time*24), sd(trip_time_new2_biting_time*24)),
                                                                          y2 = dnorm(x, mean(clusts_median_Mutasa), sd(clusts_median_Mutasa))))

clust_times_part2 <- as.data.frame(cbind(c(Mutasa5_clust_trips_shorter2_no_home_biting_places2$trip_time_new2_biting_time*24, clusts_median_Mutasa),
                                         c(rep("overall",nrow(Mutasa5_clust_trips_shorter2_no_home_biting_places2)),
                                           rep("part_median",length(clusts_median_Mutasa)))))
colnames(clust_times_part2) <- c("values","type")
clust_times_part2$values <- as.numeric(as.character(clust_times_part2$values))
clust_times_part2$type <- factor(clust_times_part2$type)

ggplot(data=clust_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,275)) +
  scale_x_continuous(breaks=c(seq(0,300,50)))

## just medians
ggplot(data=clust_times_part2[which(clust_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 5, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=2) +
  coord_cartesian(xlim=c(0,250)) +
  scale_x_continuous(breaks=c(seq(0,300,50)))
#
#
##### percent biting time spent at home vs elsewhere #####
Mutasa5_clusts_shorter2_biting_places <- Mutasa5_clusts_shorter2[which(!(is.na(Mutasa5_clusts_shorter2$clust_trips_time_new_biting_time))),]

Mutasa5_clusts_shorter2_biting_places$percent_home_biting_time <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_biting_places$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  percent <- Mutasa5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Mutasa5_clusts_shorter2_biting_places$partid == part & Mutasa5_clusts_shorter2_biting_places$home_clust2 == 1)]/
    sum(Mutasa5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Mutasa5_clusts_shorter2_biting_places$partid == part)])
  if(length(percent) == 0){
    percent <- 0
  }
  Mutasa5_clusts_shorter2_biting_places$percent_home_biting_time[which(Mutasa5_clusts_shorter2_biting_places$partid == part)] <- percent*100
}

Mutasa5_clusts_shorter2_biting_places2 <- Mutasa5_clusts_shorter2_biting_places[!duplicated(Mutasa5_clusts_shorter2_biting_places$partid),]
## none with 0%
homes <- Mutasa5_clusts_shorter2_home$partid
Mutasa5_clusts_shorter2_biting_places3 <- Mutasa5_clusts_shorter2_biting_places[which(Mutasa5_clusts_shorter2_biting_places$partid %in% homes),]
Mutasa5_clusts_shorter2_biting_places3b <- Mutasa5_clusts_shorter2_biting_places3[!duplicated(Mutasa5_clusts_shorter2_biting_places3$partid),]

favstats(Mutasa5_clusts_shorter2_biting_places3b$percent_home_biting_time)
favstats(Mutasa5_clusts_shorter2_biting_places3b$percent_home_biting_time[which(Mutasa5_clusts_shorter2_biting_places3b$percent_home_biting_time != 0)])

ggplot(data=Mutasa5_clusts_shorter2_biting_places3b, aes(x=(percent_home_biting_time))) +
  geom_histogram(binwidth = 10, fill = "black") +
  coord_cartesian(xlim=c(0,100)) +
  scale_x_continuous(breaks=c(seq(0,100,10)))
##
##
###############
###############
##### NCHELENGE #####
setwd("/Users/kschaber/Desktop/ICEMR_GPSLoggers/Data")
lake_bound <- st_read("Lake_Mweru2.shp")
lake_bound2 <- as_Spatial(lake_bound)
point <- data.frame(lon=Nchelenge5_clusts_shorter2$log_long, lat=Nchelenge5_clusts_shorter2$log_lat)
sp2   <- SpatialPoints(point,proj4string=CRS(proj4string(lake_bound2)))
Nchelenge5_clusts_shorter2$in_lake <- as.factor(as.integer(gContains(lake_bound2,sp2, byid=TRUE)))
summary(Nchelenge5_clusts_shorter2$in_lake)

zambia <- st_read("Zambia.shp")
zambia2 <- as_Spatial(zambia)
point <- data.frame(lon=Nchelenge5_clusts_shorter2$log_long, lat=Nchelenge5_clusts_shorter2$log_lat)
sp2   <- SpatialPoints(point,proj4string=CRS(proj4string(zambia2)))
Nchelenge5_clusts_shorter2$in_zambia <- as.factor(as.integer(gContains(zambia2,sp2, byid=TRUE)))
summary(Nchelenge5_clusts_shorter2$in_zambia) ## no one leaves country

setwd("/Users/kschaber/Desktop/ICEMR_GPSLoggers/Code/HDBSCAN/Outputs/New/3.17")
#######
Nchelenge5_clusts_shorter2$lakeside_home <- factor(NA, levels=c(0,1))
parts <- levels(as.factor(Nchelenge5_clusts_shorter2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  if(length(which(Nchelenge5_clusts_shorter2$partid == part & Nchelenge5_clusts_shorter2$home_clust2 == 1)) > 0){
    Nchelenge5_clusts_shorter2$lakeside_home[which(Nchelenge5_clusts_shorter2$partid == part)] <- 
      Nchelenge5_clusts_shorter2$lakeside[which(Nchelenge5_clusts_shorter2$partid == part & Nchelenge5_clusts_shorter2$home_clust2 == 1)]
  }
  else{
    Nchelenge5_clusts_shorter2$lakeside_home[which(Nchelenge5_clusts_shorter2$partid == part)] <- NA
  }
}
#### Get values for biting time period ####
Nchelenge5_new2 <- Nchelenge5_new
Nchelenge5_new2$biting_time <- factor(0, levels = c(0,1))
Nchelenge5_new2$biting_time[which(Nchelenge5_new2$time_new_hour %in% c(21:23,0:5))] <- 1

Nchelenge5_new2$trip_time_new2_biting_time <- as.numeric(NA)
Nchelenge5_new2$clust_trips_time_new_biting_time <- as.numeric(NA)

parts <- levels(as.factor(Nchelenge5_new2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Nchelenge5_new2$clust[which(Nchelenge5_new2$partid == part)]))
  for(jj in 1:length(clusts)){
    clust <- clusts[jj]
    trips <- levels(as.factor(Nchelenge5_new2$trip_count_new[which(Nchelenge5_new2$partid == part & Nchelenge5_new2$clust == clust)]))
    temp3 <- Nchelenge5_new2[which(Nchelenge5_new2$partid == part & Nchelenge5_new2$clust == clust & Nchelenge5_new2$biting_time == 1),]
    if(nrow(temp3) > 0){
      xx2 <- aggregate(time_atpoint_new ~ trip_count_new, temp3, sum)
      for(kk in 1:nrow(xx2)){
        Nchelenge5_new2$trip_time_new2_biting_time[which(Nchelenge5_new2$partid == part & 
                                                        Nchelenge5_new2$clust == clust &
                                                        Nchelenge5_new2$trip_count_new == xx2$trip_count_new[kk])] <- xx2$time_atpoint_new[kk]
      }
    }
  }
  temp4 <- Nchelenge5_new2[which(Nchelenge5_new2$partid == part),]
  temp4b <- temp4[!duplicated(temp4[c(1,21,25)]),]
  if(nrow(temp4b[which(!(is.na(temp4b$trip_time_new2_biting_time))),])>0){
    biting <- aggregate(trip_time_new2_biting_time ~ clust, temp4b, sum)
    if(nrow(biting) > 0){
      for(kk in 1:nrow(biting)){
        Nchelenge5_new2$clust_trips_time_new_biting_time[which(Nchelenge5_new2$partid == part & Nchelenge5_new2$clust == biting$clust[kk] & !(is.na(Nchelenge5_new2$trip_count_new)))] <- biting$trip_time_new2_biting_time[kk]
      }
    }
  }
}
Nchelenge5_new2$part_trips_time_new_biting_time <- as.numeric(0)
temp <- Nchelenge5_new2[!duplicated(Nchelenge5_new2[c(1,21)]),]
tot <- aggregate(clust_trips_time_new_biting_time ~ partid, temp, sum)
for(ii in 1:nrow(tot)){
  Nchelenge5_new2$part_trips_time_new_biting_time[which(Nchelenge5_new2$partid == tot$partid[ii])] <- tot$clust_trips_time_new_biting_time[ii]
}

prop.table(table(Nchelenge5_new2$home_clust2, Nchelenge5_new2$biting_time, deparse.level = 2), margin=2)*100
#                           Nchelenge5_new2$biting_time
# Nchelenge5_new2$home_clust2           0           1
# 0                           40.53787206 30.08914649
# 1                           59.46212794 69.91085351
prop.table(table(Nchelenge5_new2$home_clust2, Nchelenge5_new2$biting_time, deparse.level = 2), margin=1)*100
#                           Nchelenge5_new2$biting_time
# Nchelenge5_new2$home_clust2           0           1
# 0                           81.39869043 18.60130957
# 1                           73.42266266 26.57733734

#### 70% of observations taken during biting times were in the home cluster
#### 27% of obervations in the home cluster were taken during biting times

Nchelenge5_clusts_shorter_temp <- Nchelenge5_new2[order(Nchelenge5_new2$partid, Nchelenge5_new2$clust),]
Nchelenge5_clusts_shorter_temp2 <- Nchelenge5_clusts_shorter_temp[!duplicated(Nchelenge5_clusts_shorter_temp[,c(1,21)]),c(1:6,15,16,11,12,18:24,28,29:34,45:49,40:44,50,51:65,69,68)]

Nchelenge5_clust_trips_shorter_temp <- Nchelenge5_new2[order(Nchelenge5_new2$partid, Nchelenge5_new2$clust, Nchelenge5_new2$trip_count),]
Nchelenge5_clust_trips_shorter_temp2 <- Nchelenge5_clust_trips_shorter_temp[!duplicated(Nchelenge5_clust_trips_shorter_temp[,c(1,21,25)]),c(11:6,15,16,11,12,18:22,23,25:27,28,29:34,45:49,40:44,35:39,50,51:65,69,68,67)]


Nchelenge5_clusts_shorter2$part_trips_time_new_biting_time <- Nchelenge5_clusts_shorter_temp2$part_trips_time_new_biting_time
Nchelenge5_clusts_shorter2$clust_trips_time_new_biting_time <- Nchelenge5_clusts_shorter_temp2$clust_trips_time_new_biting_time

Nchelenge5_clust_trips_shorter2$part_trips_time_new_biting_time <- Nchelenge5_clust_trips_shorter_temp2$part_trips_time_new_biting_time
Nchelenge5_clust_trips_shorter2$clust_trips_time_new_biting_time <- Nchelenge5_clust_trips_shorter_temp2$clust_trips_time_new_biting_time
Nchelenge5_clust_trips_shorter2$trip_time_new2_biting_time <- Nchelenge5_clust_trips_shorter_temp2$trip_time_new2_biting_time
##
### number of locations (without home) ####
Nchelenge5_clusts_shorter2_no_home <- Nchelenge5_clusts_shorter2[which(Nchelenge5_clusts_shorter2$home_clust2 == 0),]
loc_counts <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)$n
only_home <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home$partid))
loc_counts2 <- c(loc_counts,rep(0, only_home))
favstats(loc_counts2)
# min Q1 median Q3 max     mean       sd  n missing
#    0 3.5      8 12  23 8.573333 5.688237 75       0
hist(loc_counts2, xlim=c(0,25))

### distance of locations (without home) ####
### TOTAL ###
Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine <- Nchelenge5_clusts_shorter2_no_home$hhdist_m_haversine2/1000
favstats(Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#  0 0.3695943 1.240013 4.580401 187.2311 6.848572 20.13687 643       0

### PER PART ###
medians_km <- (favstats(Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home$partid)$median)
favstats(medians_km) # km
#          min        Q1   median       Q3      max     mean       sd  n missing
#  0.1012507 0.5565415 1.166392 3.160273 186.7748 6.633795 24.30083 71       0

## REMOVE OUTLIER ##
medians_km_short <- medians_km[-which(medians_km > 70)]
favstats(medians_km_short)
# min        Q1    median       Q3      max     mean       sd  n missing
# 0.1012507 0.5552885 1.15471 3.006862 25.71711 2.881181 5.010192 69       0

x <- seq(0,500, length.out=500)
df <- with(Nchelenge5_clusts_shorter2_no_home, data.frame(x = x, y1 = dnorm(x, mean(hhdist_km_haversine), sd(hhdist_km_haversine)),
                                                            y2 = dnorm(x, mean(medians_km_short), sd(medians_km_short))))

dists_km2 <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine, medians_km_short),
                                 c(rep("overall",nrow(Nchelenge5_clusts_shorter2_no_home)),
                                   rep("part_median",length(medians_km_short)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,200))

## just medians
ggplot(data=dists_km2[which(dists_km2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 2, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=5) +
  coord_cartesian(xlim=c(0,30))
#
### number of trips (without home) ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2$lakeside_home <- factor(NA, levels=c(0,1))
parts <- levels(as.factor(Nchelenge5_clust_trips_shorter2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  if(length(which(Nchelenge5_clust_trips_shorter2$partid == part & Nchelenge5_clust_trips_shorter2$home_clust2 == 1)) > 0){
    Nchelenge5_clust_trips_shorter2$lakeside_home[which(Nchelenge5_clust_trips_shorter2$partid == part)] <- 
      Nchelenge5_clust_trips_shorter2$lakeside[which(Nchelenge5_clust_trips_shorter2$partid == part & Nchelenge5_clust_trips_shorter2$home_clust2 == 1)][1]
  }
  else{
    Nchelenge5_clust_trips_shorter2$lakeside_home[which(Nchelenge5_clust_trips_shorter2$partid == part)] <- NA
  }
}

Nchelenge5_clust_trips_shorter2_no_home <- Nchelenge5_clust_trips_shorter2[which(Nchelenge5_clust_trips_shorter2$home_clust2 == 0),]
Nchelenge5_clust_trips_shorter2_no_home2 <- Nchelenge5_clust_trips_shorter2_no_home[order(Nchelenge5_clust_trips_shorter2_no_home$partid,
                                                                                          Nchelenge5_clust_trips_shorter2_no_home$clust,
                                                                                          -Nchelenge5_clust_trips_shorter2_no_home$trip_count_new),]
Nchelenge5_clust_trips_shorter2_no_home3 <- Nchelenge5_clust_trips_shorter2_no_home2[!duplicated(Nchelenge5_clust_trips_shorter2_no_home2[,c(1,14)]),]

favstats(Nchelenge5_clust_trips_shorter2_no_home3$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      2  3  23 2.515528 2.712463 644       0

### PER PART ###
trips_part_median <- favstats(Nchelenge5_clust_trips_shorter2_no_home3$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3$partid)$median
favstats(trips_part_median)
#   min Q1 median Q3 max     mean       sd  n missing
#   1  1      2  2   4 1.640845 0.6163435 71       0

### PLOT BOTH ###
x <- seq(0,500, length.out=500)
df <- with(Nchelenge5_clust_trips_shorter2_no_home3, data.frame(x = x, y1 = dnorm(x, mean(trip_count_new), sd(trip_count_new)),
                                                          y2 = dnorm(x, mean(trips_part_median), sd(trips_part_median))))

trips_part2 <- as.data.frame(cbind(c(Nchelenge5_clust_trips_shorter2_no_home3$trip_count_new, trips_part_median),
                                   c(rep("overall",nrow(Nchelenge5_clust_trips_shorter2_no_home3)),
                                     rep("part_median",length(trips_part_median)))))
colnames(trips_part2) <- c("values","type")
trips_part2$values <- as.numeric(as.character(trips_part2$values))
trips_part2$type <- factor(trips_part2$type)

ggplot(data=trips_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,35))

## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=2) +
  coord_cartesian(xlim=c(0,10))+ 
  scale_x_continuous(breaks=c(seq(0,14,1)))
#
#
### time per trip ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour <- Nchelenge5_clust_trips_shorter2_no_home$trip_time_new2*24
favstats(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
#  0.2416667 0.5752778 1.182639 2.709722 394.7653 5.679629 24.06426 1610       0

### PER PART ###
trip_time_part_median <- favstats(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home$partid)$median
favstats(trip_time_part_median)
#       min        Q1   median       Q3      max     mean       sd  n missing
#  0.4844444 0.8788889 1.312222 1.854306 6.663194 1.579349 1.109373 71       0

### PLOT BOTH ###
x <- seq(0,700, length.out=5000)
df <- with(Nchelenge5_clust_trips_shorter2_no_home, data.frame(x = x, y1 = dnorm(x, mean(trip_time_new_hour), sd(trip_time_new_hour)),
                                                               y2 = dnorm(x, mean(trip_time_part_median), sd(trip_time_part_median))))

trip_times_part2 <- as.data.frame(cbind(c(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour, trip_time_part_median),
                                        c(rep("overall",nrow(Nchelenge5_clust_trips_shorter2_no_home)),
                                          rep("part_median",length(trip_time_part_median)))))
colnames(trip_times_part2) <- c("values","type")
trip_times_part2$values <- as.numeric(as.character(trip_times_part2$values))
trip_times_part2$type <- factor(trip_times_part2$type)

ggplot(data=trip_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "darkgrey") +
  coord_cartesian(xlim=c(0,400))
#
## just medians
ggplot(data=trip_times_part2[which(trip_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 1, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1.5) +
  coord_cartesian(xlim=c(0,10))+ 
  scale_x_continuous(breaks=c(seq(0,14,1)))
#
### Home Loc ####
Nchelenge5_clusts_shorter2_home <- Nchelenge5_clusts_shorter2[which(Nchelenge5_clusts_shorter2$home_clust2 == 1),]
no_home <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_home$partid))
### number of 'trips' to home loc ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2_home <- Nchelenge5_clust_trips_shorter2[which(Nchelenge5_clust_trips_shorter2$home_clust2 == 1),]
Nchelenge5_clust_trips_shorter2_home2 <- Nchelenge5_clust_trips_shorter2_home[order(Nchelenge5_clust_trips_shorter2_home$partid,
                                                                                    Nchelenge5_clust_trips_shorter2_home$clust,
                                                                                    -Nchelenge5_clust_trips_shorter2_home$trip_count_new),]
Nchelenge5_clust_trips_shorter2_home3 <- Nchelenge5_clust_trips_shorter2_home2[!duplicated(Nchelenge5_clust_trips_shorter2_home2[,c(1,14)]),]

favstats(Nchelenge5_clust_trips_shorter2_home3$trip_count_new)
#   min   Q1 median    Q3  max  mean       sd  n missing
# 1  5    9.5 19.75  60 14.68919 13.90272 74       0

### total time at home loc ####
Nchelenge5_clust_trips_shorter2_home4 <- Nchelenge5_clust_trips_shorter2_home[!duplicated(Nchelenge5_clust_trips_shorter2_home$partid),]
favstats(Nchelenge5_clust_trips_shorter2_home4$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
#   0.4103125 17.89881 24.92939 29.88956 40.7908 23.55836 8.936888 74       0
#
### percent time at home ####
Nchelenge5_clusts_shorter2_home$perc_time_home <- (Nchelenge5_clusts_shorter2_home$clust_trips_time_new/Nchelenge5_clusts_shorter2_home$part_trips_time_new)*100

favstats(Nchelenge5_clusts_shorter2_home$perc_time_home)
#  min       Q1   median       Q3   max    mean       sd  n missing
# 1.618541 81.32746 91.74501 98.0249 100 84.15788 20.95541 74       0

Nchelenge5_clusts_shorter2_home$partid[which(Nchelenge5_clusts_shorter2_home$perc_time_home < 5)]
## "22030"

ggplot(data=Nchelenge5_clusts_shorter2_home, aes(x=(perc_time_home))) +
  geom_histogram(binwidth = 10, fill = "black") +
  coord_cartesian(xlim=c(0,100))
#
# examine those with low percent time at home #####
parts_low_perc_home <- Nchelenge5_clusts_shorter2_home$partid[which(Nchelenge5_clusts_shorter2_home$perc_time_home < 15)]
temp_trips <- Nchelenge5_clust_trips_shorter2[which(Nchelenge5_clust_trips_shorter2$partid %in% parts_low_perc_home),]
temp_clusts <- Nchelenge5_clusts_shorter2[which(Nchelenge5_clusts_shorter2$partid %in% parts_low_perc_home),]
temp_clusts_no_home <- temp_clusts[which(temp_clusts$home_clust2 == 0),]

### number of locs
favstats(favstats(temp_clusts_no_home$clust~temp_clusts_no_home$partid)$n)
# min Q1  median    Q3 max     mean    sd    n  missing
#   9 10.5     12 12.5  13 11.33333 2.081666 3       0
# VS overall:
# min Q1 median Q3 max     mean       sd  n missing
#   1  7     12 15  17 10.77419 4.267373 62       0
hist(loc_counts2, xlim=c(0,25))
hist(favstats(temp_clusts_no_home$clust~temp_clusts_no_home$partid)$n, xlim=c(0,25), add=TRUE, col="red")
#
#### dist of locs
medians_km_all <- (favstats(Nchelenge5_clusts_shorter2_no_home$hhdist_m_haversine~Nchelenge5_clusts_shorter2_no_home$partid)$median)/1000
medians_km_low <- (favstats(temp_clusts_no_home$hhdist_m_haversine~temp_clusts_no_home$partid)$median)/1000
favstats(medians_km_low)
#    min    Q1 median    Q3   max  mean    sd  n missing
# 1.312188 2.413546 3.514903 4.46348 5.412056 3.413049 2.051831 3       0

x <- seq(0,500, length.out=5000)
df <- data.frame(x = x, y1 = dnorm(x, mean(medians_km_all), sd(medians_km_all)),
                 y2 = dnorm(x, mean(medians_km_low), sd(medians_km_low)))
dists_km2 <- as.data.frame(cbind(c(medians_km_all, medians_km_low),
                                 c(rep("all",length(medians_km_all)),
                                   rep("low perc",length(medians_km_low)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) + 
  scale_fill_manual(values = c("grey","red")) +
  scale_color_manual(values = c("grey","red")) +
  coord_cartesian(xlim=c(0,30)) +
  scale_x_continuous(breaks=c(seq(0,30,5)))
##
### number of trips to locs ##
temp_trips_no_home3 <- Nchelenge5_clust_trips_shorter2_no_home3[which(Nchelenge5_clust_trips_shorter2_no_home3$partid %in% parts_low_perc_home),]
trips_part_median_all <- favstats(Nchelenge5_clust_trips_shorter2_no_home3$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3$partid)$median
trips_part_median_low <- favstats(temp_trips_no_home3$trip_count_new ~ temp_trips_no_home3$partid)$median
favstats(trips_part_median_low)
# min  Q1 median Q3 max  mean    sd  n missing
#   1  1      1 1.5   2 1.333333 0.5773503 3       0

x <- seq(1,500, length.out=5000)
df <- data.frame(x = x, y1 = dnorm(x, mean(trips_part_median_all), sd(trips_part_median_all)),
                 y2 = dnorm(x, mean(trips_part_median_low), sd(trips_part_median_low)))
dists_km2 <- as.data.frame(cbind(c(trips_part_median_all, trips_part_median_low),
                                 c(rep("all",length(trips_part_median_all)),
                                   rep("low_perc",length(trips_part_median_low)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("grey","red")) +
  scale_color_manual(values = c("grey","red")) +
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks=c(seq(0,10,1)))
##
### time per trip ##
temp_trips_no_home <- Nchelenge5_clust_trips_shorter2_no_home[which(Nchelenge5_clust_trips_shorter2_no_home$partid %in% parts_low_perc_home),]
trip_time_part_median_all <- favstats(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home$partid)$median
trip_time_part_median_low <- favstats(temp_trips_no_home$trip_time_new_hour ~ temp_trips_no_home$partid)$median
favstats(trip_time_part_median_low)
# min    Q1 median    Q3   max  mean    sd  n missing
# 1.68625 1.709653 1.733056 1.938819 2.144583 1.85463 0.2521955 3       0

x <- seq(0,500, length.out=5000)
df <- data.frame(x = x, y1 = dnorm(x, mean(trip_time_part_median_all), sd(trip_time_part_median_all)),
                 y2 = dnorm(x, mean(trip_time_part_median_low), sd(trip_time_part_median_low)))
dists_km2 <- as.data.frame(cbind(c(trip_time_part_median_all, trip_time_part_median_low),
                                 c(rep("all",length(trip_time_part_median_all)),
                                   rep("low_perc",length(trip_time_part_median_low)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) +
  scale_fill_manual(values = c("grey","red")) +
  scale_color_manual(values = c("grey","red")) +
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks=c(seq(0,10,1)))
##
### Average Direction of Travel #####
Nchelenge5_clusts_shorter2_no_home$dir_from_home_bearing <- numeric(length=length(Nchelenge5_clusts_shorter2_no_home$partid))

parts <- levels(as.factor(Nchelenge5_clusts_shorter2_no_home$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Nchelenge5_clusts_shorter2_no_home[which(Nchelenge5_clusts_shorter2_no_home$partid == part),]
  for(j in 1:nrow(temp)){
    temp$dir_from_home_deg[j] <- bearing(c(temp$hh_long[j], temp$hh_lat[j]), c(temp$log_long[j], temp$log_lat[j]))
    if(temp$dir_from_home_deg[j] < 0){
      temp$dir_from_home_deg[j] <- 360+temp$dir_from_home_deg[j]
    }
  }
  Nchelenge5_clusts_shorter2_no_home$dir_from_home_bearing[which(Nchelenge5_clusts_shorter2_no_home$partid == part)] <- temp$dir_from_home_deg
}
## to get point on unit circle, must have angle theta from
## positive x axis, going positive to the counterclockwise direction
## these are in degrees from positive y going clockwise (bearings)
Nchelenge5_clusts_shorter2_no_home$dir_from_home_deg <- 90 - Nchelenge5_clusts_shorter2_no_home$dir_from_home_bearing
Nchelenge5_clusts_shorter2_no_home$dir_from_home_deg[which(Nchelenge5_clusts_shorter2_no_home$dir_from_home_deg < 0)] <- Nchelenge5_clusts_shorter2_no_home$dir_from_home_deg[which(Nchelenge5_clusts_shorter2_no_home$dir_from_home_deg < 0)] + 360
Nchelenge5_clusts_shorter2_no_home$dir_from_home_rad <- deg2rad(Nchelenge5_clusts_shorter2_no_home$dir_from_home_deg)

## average degree and magnitude
Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_deg <- as.numeric(NA)
Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_no_home$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Nchelenge5_clusts_shorter2_no_home[which(Nchelenge5_clusts_shorter2_no_home$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_deg[which(Nchelenge5_clusts_shorter2_no_home$partid == part)] <- avg_deg
  Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_mag[which(Nchelenge5_clusts_shorter2_no_home$partid == part)] <- avg_vec_dist/1000
}

Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_bearing <- 90 - Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_deg
Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_bearing[which(Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_bearing < 0)] <- Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_bearing[which(Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_bearing < 0)] + 360


temp <- Nchelenge5_clusts_shorter2_no_home[!duplicated(Nchelenge5_clusts_shorter2_no_home$partid),]

ggplot(data=temp, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(fill="black", col="black", position=position_dodge2(preserve = "single"), width=1) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw()
#

summary(temp$avg_dir_from_home_bearing)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3489 130.8408 192.8555 200.4319 289.8498 358.9661 
summary(temp$avg_dir_from_home_mag)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1035   0.6148   1.5248   5.8709   3.5179 105.1509 

temp2 <- temp[which(temp$avg_dir_from_home_mag < 65),] ## 69 of 71 parts

summary(temp2$avg_dir_from_home_bearing)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3489 130.4685 208.9877 201.6504 290.2358 358.9661 

summary(temp2$avg_dir_from_home_mag)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1035  0.6055  1.4661  3.5150  3.2342 21.5691 

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(fill="black", col="black", position=position_dodge2(preserve = "single"), width=1) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw()




## average degree only, not accounting for distance from home
Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_no_home$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Nchelenge5_clusts_shorter2_no_home[which(Nchelenge5_clusts_shorter2_no_home$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_deg2[which(Nchelenge5_clusts_shorter2_no_home$partid == part)] <- avg_deg
}

Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_values <- 90 - Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_deg2
Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_values[which(Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_values < 0)] <- Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_values[which(Nchelenge5_clusts_shorter2_no_home$avg_dir_from_home_values < 0)] + 360
temp <- Nchelenge5_clusts_shorter2_no_home[!duplicated(Nchelenge5_clusts_shorter2_no_home$partid),]

temp$avg_dir_from_home_values[which(temp$avg_dir_from_home_values > 348.75)] <- -(360-temp$avg_dir_from_home_values[which(temp$avg_dir_from_home_values > 348.75)])
ggplot(data=temp, aes(x=avg_dir_from_home_values, y=..count..)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw()
#
### time spent at clusters vs elsewhere (travel, less important places, etc) ####
Nchelenge5_clusts_shorter2_temp <- Nchelenge5_clusts_shorter2[!duplicated(Nchelenge5_clusts_shorter2$partid),]
Nchelenge5_clusts_shorter2_temp <- Nchelenge5_clusts_shorter2_temp[which(Nchelenge5_clusts_shorter2_temp$partid != "30240"),]
Nchelenge5_clusts_shorter2_temp$perc_time_in_all_clusts <- Nchelenge5_clusts_shorter2_temp$part_trips_time_new/Nchelenge5_clusts_shorter2_temp$total_time

# percent time in clusters
summary(Nchelenge5_clusts_shorter2_temp$perc_time_in_all_clusts)
#    Min. 1st Qu.  Median    Mean  3rd Qu.    Max.
# 0.7237  0.9337  0.9670  0.9515  0.9904  1.3064 

## only issue is 22826 because added about 13% worth of missing gps hours, but spread throughout

# time outside clusters
summary(Nchelenge5_clusts_shorter2_temp$total_time - Nchelenge5_clusts_shorter2_temp$part_trips_time_new)
#       Min.  1st Qu.  Median      Mean    3rd Qu.  Max. 
#     -6.7406  0.2351  0.8802  1.4571  2.0584  9.1276 
#
### Time at locations ####
favstats(Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new*24)
# min       Q1   median       Q3      max     mean       sd   n missing
# 0.2425 1.050556 2.247778 5.035694 784.0228 14.19519 58.03024 643       0

clusts_median_Nchelenge <- favstats((Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new*24)~Nchelenge5_clusts_shorter2_no_home$partid)$median
favstats(clusts_median_Nchelenge)
# min     Q1   median       Q3      max     mean       sd  n missing
# 0.4605556 1.5625 2.337639 3.213472 17.23194 2.969454 2.708078 71       0

clust_times_part2 <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new*24, clusts_median_Nchelenge),
                                         c(rep("overall",nrow(Nchelenge5_clusts_shorter2_no_home)),
                                           rep("part_median",length(clusts_median_Nchelenge)))))
colnames(clust_times_part2) <- c("values","type")
clust_times_part2$values <- as.numeric(as.character(clust_times_part2$values))
clust_times_part2$type <- factor(clust_times_part2$type)

ggplot(data=clust_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  #geom_line(data=df, aes(x=x, y=y1), color = "black") +
  #geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,800))

## just medians
ggplot(data=clust_times_part2[which(clust_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 1, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,20)) +
  scale_x_continuous(breaks=c(seq(0,20,2)))
#
#
####### Biting time metrics #####
##### at home ####
length(unique(Nchelenge5_clusts_shorter2$partid)) # 75
length(unique(Nchelenge5_clusts_shorter2$partid[which(!(is.na(Nchelenge5_clusts_shorter2$clust_trips_time_new_biting_time)))])) # 75
## all parts spend time somewhere during biting hours
nrow(Nchelenge5_clusts_shorter2_home) # 74
Nchelenge5_clusts_shorter2_home_biting_places <- Nchelenge5_clusts_shorter2_home[which(!(is.na(Nchelenge5_clusts_shorter2_home$clust_trips_time_new_biting_time))),]
nrow(Nchelenge5_clusts_shorter2_home_biting_places) # 73
## of 74 with home location, all but 1 spend time there during biting hours

favstats(Nchelenge5_clusts_shorter2_home_biting_places$clust_trips_time_new)
# min          Q1 median          Q3         max        mean          sd  n missing
# 0.4103 17.91  24.94 29.96 40.79 23.82 8.705 73       0
favstats(Nchelenge5_clusts_shorter2_home_biting_places$clust_trips_time_new_biting_time)
# min          Q1      median          Q3         max        mean          sd  n missing
# 0.01547 6.479  7.911 10.18 14.5 8.026 3.107 73       0

# percent of time at home that is during biting hours
favstats((Nchelenge5_clusts_shorter2_home_biting_places$clust_trips_time_new_biting_time/Nchelenge5_clusts_shorter2_home_biting_places$clust_trips_time_new)*100)
# min          Q1      median          Q3         max        mean          sd  n missing
# 3.77 31.02  34.16 37.69 52.87 33.71 7.386 73       0
##
##### outside home ####
### number locations ####
Nchelenge5_clusts_shorter2_no_home_biting_places <- Nchelenge5_clusts_shorter2_no_home[which(!(is.na(Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new_biting_time))),]
## only when >30 minutes during biting time
Nchelenge5_clusts_shorter2_no_home_biting_places2 <- Nchelenge5_clusts_shorter2_no_home_biting_places[which((Nchelenge5_clusts_shorter2_no_home_biting_places$clust_trips_time_new_biting_time*24) >= 0.5),]
## removes 30 places, 3 participants

loc_counts <- favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2$partid)$n
only_home <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2$partid))
loc_counts2 <- c(loc_counts,rep(0, only_home))
favstats(loc_counts2)
# min Q1 median Q3 max     mean       sd  n missing
#    0  0      1  2  12 1.72 2.299 75       0
summary(as.factor(loc_counts2))
#  0  1  2  3  4  5  6  8  9 12 
# 28 18 12  5  5  1  3  1  1  1 
## all those with 0 have biting time at home location (not just nowhere)
hist(loc_counts2, xlim=c(-1,12), breaks=c(seq(-1,12,by=1)))
#
### time per location ####
favstats((Nchelenge5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5112 1.641  6.895 17.68 285.4 20.65 40.23 129       0

clusts_median_Nchelenge <- (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2$partid)$median)*24
favstats(clusts_median_Nchelenge)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.8406 2.558    8.5 23.51 82.97 16.83 19.74 47       0

x <- seq(0,1000, length.out=1000)
df <- with(Nchelenge5_clusts_shorter2_no_home_biting_places2, data.frame(x = x, y1 = dnorm(x, mean(clust_trips_time_new_biting_time*24), sd(clust_trips_time_new_biting_time*24)),
                                                                      y2 = dnorm(x, mean(clusts_median_Nchelenge), sd(clusts_median_Nchelenge))))

clust_times_part2 <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time*24, clusts_median_Nchelenge),
                                         c(rep("overall",nrow(Nchelenge5_clusts_shorter2_no_home_biting_places2)),
                                           rep("part_median",length(clusts_median_Nchelenge)))))
colnames(clust_times_part2) <- c("values","type")
clust_times_part2$values <- as.numeric(as.character(clust_times_part2$values))
clust_times_part2$type <- factor(clust_times_part2$type)

ggplot(data=clust_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,300))

## just medians
ggplot(data=clust_times_part2[which(clust_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 5, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,85)) +
  scale_x_continuous(breaks=c(seq(0,100,20)))
#
#
### number trips ####
Nchelenge5_clust_trips_shorter2_no_home_biting_places <- Nchelenge5_clust_trips_shorter2_no_home[which(!(is.na(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new2_biting_time))),]
## only when >30 minutes during biting time
Nchelenge5_clust_trips_shorter2_no_home_biting_places2 <- Nchelenge5_clust_trips_shorter2_no_home_biting_places[which((Nchelenge5_clust_trips_shorter2_no_home_biting_places$trip_time_new2_biting_time*24) >= 0.5),]
## removes 74 trips, 6 participants

Nchelenge5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2 <- numeric(length=nrow(Nchelenge5_clust_trips_shorter2_no_home_biting_places2))
parts <- levels(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$clust[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$partid == part)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    counts <- Nchelenge5_clust_trips_shorter2_no_home_biting_places2$trip_count_new[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$partid == part & Nchelenge5_clust_trips_shorter2_no_home_biting_places2$clust == clust)]
    counts <- as.factor(counts)
    levels(counts) <- c(1:nlevels(counts))
    Nchelenge5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$partid == part & Nchelenge5_clust_trips_shorter2_no_home_biting_places2$clust == clust)] <- as.numeric(as.character(counts))
  }
}

Nchelenge5_clust_trips_shorter2_no_home_biting_places2b <- Nchelenge5_clust_trips_shorter2_no_home_biting_places2[order(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$partid,
                                                                                                                  Nchelenge5_clust_trips_shorter2_no_home_biting_places2$clust,
                                                                                                                  -Nchelenge5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2),]
Nchelenge5_clust_trips_shorter2_no_home_biting_places3 <- Nchelenge5_clust_trips_shorter2_no_home_biting_places2b[!duplicated(Nchelenge5_clust_trips_shorter2_no_home_biting_places2b[,c(1,15)]),]

favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
# 1  1      1  2  13 2.103 2.205 126       0
summary(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2))
# 1  2  3  4  5  6  7  9 10 12 13 
# 77 22 12  3  2  2  3  2  1  1  1 

### PER PART ###
trip_count_median <- favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3$partid)$median
favstats(trip_count_median)
# min Q1 median Q3 max  mean  sd  n missing
# 1  1      1  2   7 1.766 1.22 47       0

### PLOT BOTH ###
x <- seq(0,500, length.out=5000)
df <- with(Nchelenge5_clust_trips_shorter2_no_home_biting_places3, data.frame(x = x, y1 = dnorm(x, mean(trip_count_new2), sd(trip_count_new2)),
                                                                           y2 = dnorm(x, mean(trip_count_median), sd(trip_count_median))))

trips_part2 <- as.data.frame(cbind(c(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2, trip_count_median),
                                   c(rep("overall",nrow(Nchelenge5_clust_trips_shorter2_no_home_biting_places3)),
                                     rep("part_median",length(trip_count_median)))))
colnames(trips_part2) <- c("values","type")
trips_part2$values <- as.numeric(as.character(trips_part2$values))
trips_part2$type <- factor(trips_part2$type)

ggplot(data=trips_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,13))

## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,8))+ 
  scale_x_continuous(breaks=c(seq(0,15,2)))
#
### time per trip ####
favstats((Nchelenge5_clust_trips_shorter2_no_home_biting_places2$trip_time_new2_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5031 1.132  6.414 8.994 117.1 9.991 16.46 265       0

clusts_median_Nchelenge <- (favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$trip_time_new2_biting_time~Nchelenge5_clust_trips_shorter2_no_home_biting_places2$partid)$median)*24
favstats(clusts_median_Nchelenge)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.6003 1.706  6.568 8.575 33.89 6.451 5.919 47       0

x <- seq(0,1000, length.out=1000)
df <- with(Nchelenge5_clust_trips_shorter2_no_home_biting_places2, data.frame(x = x, y1 = dnorm(x, mean(trip_time_new2_biting_time*24), sd(trip_time_new2_biting_time*24)),
                                                                           y2 = dnorm(x, mean(clusts_median_Nchelenge), sd(clusts_median_Nchelenge))))

clust_times_part2 <- as.data.frame(cbind(c(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$trip_time_new2_biting_time*24, clusts_median_Nchelenge),
                                         c(rep("overall",nrow(Nchelenge5_clust_trips_shorter2_no_home_biting_places2)),
                                           rep("part_median",length(clusts_median_Nchelenge)))))
colnames(clust_times_part2) <- c("values","type")
clust_times_part2$values <- as.numeric(as.character(clust_times_part2$values))
clust_times_part2$type <- factor(clust_times_part2$type)

ggplot(data=clust_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  scale_fill_manual(values = c("black", "grey")) +
  geom_line(data=df, aes(x=x, y=y1), color = "black") +
  geom_line(data=df, aes(x=x, y=y2), color = "grey") +
  coord_cartesian(xlim=c(0,120)) +
  scale_x_continuous(breaks=c(seq(0,300,20)))

## just medians
ggplot(data=clust_times_part2[which(clust_times_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 2, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,35)) +
  scale_x_continuous(breaks=c(seq(0,50,5)))
#
#
##### percent biting time spent at home vs elsewhere #####
Nchelenge5_clusts_shorter2_biting_places <- Nchelenge5_clusts_shorter2[which(!(is.na(Nchelenge5_clusts_shorter2$clust_trips_time_new_biting_time))),]

Nchelenge5_clusts_shorter2_biting_places$percent_home_biting_time <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_biting_places$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  percent <- Nchelenge5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Nchelenge5_clusts_shorter2_biting_places$partid == part & Nchelenge5_clusts_shorter2_biting_places$home_clust2 == 1)]/
    sum(Nchelenge5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Nchelenge5_clusts_shorter2_biting_places$partid == part)])
  if(length(percent) == 0){
    percent <- 0
  }
  Nchelenge5_clusts_shorter2_biting_places$percent_home_biting_time[which(Nchelenge5_clusts_shorter2_biting_places$partid == part)] <- percent*100
}

Nchelenge5_clusts_shorter2_biting_places2 <- Nchelenge5_clusts_shorter2_biting_places[!duplicated(Nchelenge5_clusts_shorter2_biting_places$partid),]
## 2 with 0%, but 1 of these doesn't have home location
homes <- Nchelenge5_clusts_shorter2_home$partid
Nchelenge5_clusts_shorter2_biting_places3 <- Nchelenge5_clusts_shorter2_biting_places[which(Nchelenge5_clusts_shorter2_biting_places$partid %in% homes),]
Nchelenge5_clusts_shorter2_biting_places3b <- Nchelenge5_clusts_shorter2_biting_places3[!duplicated(Nchelenge5_clusts_shorter2_biting_places3$partid),]

favstats(Nchelenge5_clusts_shorter2_biting_places3b$percent_home_biting_time)
favstats(Nchelenge5_clusts_shorter2_biting_places3b$percent_home_biting_time[which(Nchelenge5_clusts_shorter2_biting_places3b$percent_home_biting_time != 0)])

ggplot(data=Nchelenge5_clusts_shorter2_biting_places3b, aes(x=(percent_home_biting_time))) +
  geom_histogram(binwidth = 10, fill = "black") +
  coord_cartesian(xlim=c(0,100)) +
  scale_x_continuous(breaks=c(seq(0,100,10)))
##
##
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
###### ###### GENDER SPLIT ###### ######
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
###############
###############
##### CHOMA #####
###### General Gender Analysis ####
Choma_part_info <- Choma5_clusts_shorter2[!duplicated(Choma5_clusts_shorter2$partid),c(1:3)]
summary(as.factor(Choma_part_info$male)) 
#   0   1
#  31  31 
Choma_part_info$partid2 <- gsub("b","", Choma_part_info$partid)
Choma_part_info2 <- Choma_part_info[!duplicated(Choma_part_info$partid2),]
summary(as.factor(Choma_part_info2$male))
###### Split by Gender ######
### number of locations (without home) ####
Choma5_clusts_shorter2_no_home_female <- Choma5_clusts_shorter2_no_home[which(Choma5_clusts_shorter2_no_home$male == 0),]
Choma5_clusts_shorter2_no_home_male <- Choma5_clusts_shorter2_no_home[which(Choma5_clusts_shorter2_no_home$male == 1),]
loc_counts_female <- favstats(Choma5_clusts_shorter2_no_home_female$clust~Choma5_clusts_shorter2_no_home_female$partid)$n
loc_counts_male <- favstats(Choma5_clusts_shorter2_no_home_male$clust~Choma5_clusts_shorter2_no_home_male$partid)$n
only_home_female <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_female$partid))
only_home_male <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_male$partid))
loc_counts_female2 <- c(loc_counts_female,rep(0, only_home_female))
loc_counts_male2 <- c(loc_counts_male,rep(0, only_home_male))
favstats(loc_counts_female2)
# min   Q1 median     Q3  max     mean       sd  n missing
#  0 5.25   11.5 14.5  22 10.73333 5.854166 30       0
favstats(loc_counts_male2)
# min Q1 median Q3 max     mean       sd  n missing
# 3 8.5     13 18  21 13.16129 5.489971 31       0

locs_count2 <- as.data.frame(cbind(c(loc_counts_female2, loc_counts_male2),
                                   c(rep("female",length(loc_counts_female2)),
                                     rep("male",length(loc_counts_male2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

kruskal.test(locs_count2$values, locs_count2$type) # p = 0.12

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-2,25)) +
  scale_x_continuous(breaks=c(seq(0,25,5)))
#
### distance of locations (without home) ####
### TOTAL ###
favstats(Choma5_clusts_shorter2_no_home_female$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#  0.03054886 1.236443 11.21437 16.95856 112.2062 12.17021 14.48056 322       0
favstats(Choma5_clusts_shorter2_no_home_male$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 0.06883561 1.059231 3.639848 14.94968 79.10662 9.772464 12.8633 408       0

### PER PART ###
medians_km_female <- (favstats(Choma5_clusts_shorter2_no_home_female$hhdist_m_haversine~Choma5_clusts_shorter2_no_home_female$partid)$median)/1000
medians_km_male <- (favstats(Choma5_clusts_shorter2_no_home_male$hhdist_m_haversine~Choma5_clusts_shorter2_no_home_male$partid)$median)/1000
favstats(medians_km_female) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.2448758 0.7584151 7.701106 14.78895 22.63018 8.870115 7.482457 29       0
favstats(medians_km_male) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#     0.4184579 0.9994479 2.0989 12.94367 32.16992 7.310015 8.496904 31       0

dists_km2_gender <- as.data.frame(cbind(c(medians_km_female, medians_km_male),
                                        c(rep("female",length(medians_km_female)),
                                          rep("male",length(medians_km_male)))))
colnames(dists_km2_gender) <- c("values","type")
dists_km2_gender$values <- as.numeric(as.character(dists_km2_gender$values))
dists_km2_gender$type <- factor(dists_km2_gender$type)

kruskal.test(dists_km2_gender$values, dists_km2_gender$type) # p = 0.12

ggplot(data=dists_km2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,32))
#
### number of trips (without home) ####
### TOTAL ###
Choma5_clust_trips_shorter2_no_home3_female <- Choma5_clust_trips_shorter2_no_home3[which(Choma5_clust_trips_shorter2_no_home3$male == 0),]
Choma5_clust_trips_shorter2_no_home3_male <- Choma5_clust_trips_shorter2_no_home3[which(Choma5_clust_trips_shorter2_no_home3$male == 1),]

favstats(Choma5_clust_trips_shorter2_no_home3_female$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      1  3  31 2.956522 4.305105 322       0
favstats(Choma5_clust_trips_shorter2_no_home3_male$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      2  3  58 3.178922 5.982353 408       0

### PER PART ###
trips_part_median_female <- favstats(Choma5_clust_trips_shorter2_no_home3_female$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_female$partid)$median
trips_part_median_male <- favstats(Choma5_clust_trips_shorter2_no_home3_male$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_male$partid)$median
favstats(trips_part_median_female)
#   min  Q1   median   Q3      max       mean       sd  n missing
#    1  1    1.5  2   4 1.603448 0.7243138 29       0
favstats(trips_part_median_male)
#   min Q1 median  Q3   max      mean       sd  n missing
#   1  1      1  2   3  1.5 0.6191392 31       0

trips_part2_gender <- as.data.frame(cbind(c(trips_part_median_female, trips_part_median_male),
                                          c(rep("female",length(trips_part_median_female)),
                                            rep("male",length(trips_part_median_male)))))
colnames(trips_part2_gender) <- c("values","type")
trips_part2_gender$values <- as.numeric(as.character(trips_part2_gender$values))
trips_part2_gender$type <- factor(trips_part2_gender$type)

kruskal.test(trips_part2_gender$values, trips_part2_gender$type) # p = 0.7

ggplot(data=trips_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,10)) + 
  scale_x_continuous(breaks=c(seq(0,10,2)))

#
### time per trip ####
### TOTAL ###
Choma5_clust_trips_shorter2_no_home_female <- Choma5_clust_trips_shorter2_no_home[which(Choma5_clust_trips_shorter2_no_home$male == 0),]
Choma5_clust_trips_shorter2_no_home_male <- Choma5_clust_trips_shorter2_no_home[which(Choma5_clust_trips_shorter2_no_home$male == 1),]
favstats(Choma5_clust_trips_shorter2_no_home_female$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2402778 0.7197222 1.660278 5.354722 456.055 12.0874 39.30908 940       0
favstats(Choma5_clust_trips_shorter2_no_home_male$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2402778 0.6345833 1.319306 3.241042 412.5603 6.474257 22.21797 1286       0

### PER PART ###
trip_time_part_median_female <- favstats(Choma5_clust_trips_shorter2_no_home_female$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_female$partid)$median
trip_time_part_median_male <- favstats(Choma5_clust_trips_shorter2_no_home_male$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_male$partid)$median
favstats(trip_time_part_median_female)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.6075 1.2 1.545139 2.4 5.650556 1.88568 1.133672 29       0
favstats(trip_time_part_median_male)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.5879167 1.087153 1.439722 1.523958 4.800139 1.417894 0.7264819 31       0

trip_times_part2_gender <- as.data.frame(cbind(c(trip_time_part_median_female, trip_time_part_median_male),
                                               c(rep("female",length(trip_time_part_median_female)),
                                                 rep("male",length(trip_time_part_median_male)))))
colnames(trip_times_part2_gender) <- c("values","type")
trip_times_part2_gender$values <- as.numeric(as.character(trip_times_part2_gender$values))
trip_times_part2_gender$type <- factor(trip_times_part2_gender$type)
kruskal.test(trip_times_part2_gender$values, trip_times_part2_gender$type) # p = 0.7

ggplot(data=trip_times_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 0.5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,6)) +
  scale_x_continuous(breaks=c(seq(0,6,1)))
#
#
### Home Loc ####
Choma5_clusts_shorter2_home_female <- Choma5_clusts_shorter2_home[which(Choma5_clusts_shorter2_home$male == 0),] 
Choma5_clusts_shorter2_home_male <- Choma5_clusts_shorter2_home[which(Choma5_clusts_shorter2_home$male == 1),] 

no_home_female <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Choma5_clusts_shorter2_home_female$partid))# 1
## 1 female with no home
# # "15927"
no_home_male <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Choma5_clusts_shorter2_home_male$partid))# 3
## 2 males with no homes 
# # "10692b" "10753"

### number of 'trips' to home loc ####
### TOTAL ###
Choma5_clust_trips_shorter2_home3_female <- Choma5_clust_trips_shorter2_home3[which(Choma5_clust_trips_shorter2_home3$male == 0),]
Choma5_clust_trips_shorter2_home3_male <- Choma5_clust_trips_shorter2_home3[which(Choma5_clust_trips_shorter2_home3$male == 1),]

favstats(Choma5_clust_trips_shorter2_home3_female$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#    1  4      5  9  38    9 10.00714 29       0
favstats(Choma5_clust_trips_shorter2_home3_male$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#     2  4     10 19  37 13.58621 11.18838 29       0

### total time at home loc ####
Choma5_clust_trips_shorter2_home4_female <- Choma5_clust_trips_shorter2_home4[which(Choma5_clust_trips_shorter2_home4$male == 0),]
Choma5_clust_trips_shorter2_home4_male <- Choma5_clust_trips_shorter2_home4[which(Choma5_clust_trips_shorter2_home4$male == 1),]

favstats(Choma5_clust_trips_shorter2_home4_female$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
# 0.02523149 2.456152 18.52784 28.67421 116.9807 20.53227 22.66617 29       0
favstats(Choma5_clust_trips_shorter2_home4_male$clust_trips_time_new)
#         min       Q1   median       Q3      max     mean       sd  n missing
# 0.0630382 14.46659 23.30234 24.74807 32.9477 19.24256 10.16361 29       0
#
### percent time at home ####
## without 3 people who spent 0% of time at home
favstats(Choma5_clusts_shorter2_home_female$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
#   0.08520987 7.064713 75.7409 96.04384 100 54.86957 41.52897 29       0
favstats(Choma5_clusts_shorter2_home_male$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#   0.243299 32.15271 85.37115 94.76292 99.27781 66.60755 35.88741 29       0

perc_time_home2_gender <- as.data.frame(cbind(c(Choma5_clusts_shorter2_home_female$perc_time_home, Choma5_clusts_shorter2_home_male$perc_time_home),
                                              c(rep("female",length(Choma5_clusts_shorter2_home_female$perc_time_home)),
                                                rep("male",length(Choma5_clusts_shorter2_home_male$perc_time_home)))))
colnames(perc_time_home2_gender) <- c("values","type")
perc_time_home2_gender$values <- as.numeric(as.character(perc_time_home2_gender$values))
perc_time_home2_gender$type <- factor(perc_time_home2_gender$type)

kruskal.test(perc_time_home2_gender$values, perc_time_home2_gender$type) # p = 0.26

ggplot(data=perc_time_home2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,100))
##
### Average Direction of Travel #####
Choma5_clusts_shorter2_no_home$male <- as.factor(Choma5_clusts_shorter2_no_home$male)
temp <- Choma5_clusts_shorter2_no_home[!duplicated(Choma5_clusts_shorter2_no_home$partid),]
temp2 <- temp[which(!(is.na(temp$male))),]

favstats(temp2$avg_dir_from_home_bearing~temp2$male)
# temp2$male      min       Q1    median       Q3      max     mean        sd  n missing
# 1          0  8.43748 53.88350  88.09797 181.0919 305.7516 129.4633  97.04688 29       0
# 2          1 21.03243 75.74228 141.73065 272.5941 344.5980 163.6677 104.28661 31       0
favstats(temp2$avg_dir_from_home_mag~temp2$male)
# temp2$male       min        Q1    median       Q3      max     mean       sd  n missing
# 1          0 0.3610027 1.3928157 11.659669 14.83266 33.03357 9.594053 8.540664 29       0
# 2          1 0.1915475 0.9597069  3.421198 13.90912 30.09577 7.952992 8.326652 31       0

kruskal.test(temp2$avg_dir_from_home_bearing, temp2$male)

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=male), position=position_dodge2(preserve = "single"), width=2) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Male", values=c("red","blue"))




temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)] <- -(360-temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)])
ggplot(data=temp2, aes(x=avg_dir_from_home_values, y=..count.., fill=male)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~male, labeller = labeller(male = c("0"="Female", "1"="Male")))+
  scale_fill_manual(values=c("red","blue"))

#
#

### time spent at clusters vs elsewhere (travel, less important places, etc) ####
Choma5_clusts_shorter2_female <- Choma5_clusts_shorter2[which(Choma5_clusts_shorter2$male == 0),]
Choma5_clusts_shorter2_male <- Choma5_clusts_shorter2[which(Choma5_clusts_shorter2$male == 1),]
Choma5_clusts_shorter2_female$perc_time_in_all_clusts <- Choma5_clusts_shorter2_female$part_trips_time_new/Choma5_clusts_shorter2_female$total_time
Choma5_clusts_shorter2_male$perc_time_in_all_clusts <- Choma5_clusts_shorter2_male$part_trips_time_new/Choma5_clusts_shorter2_male$total_time

summary(Choma5_clusts_shorter2_female$perc_time_in_all_clusts[!duplicated(Choma5_clusts_shorter2_female$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7556  0.8992  0.9560  0.9374  0.9872  0.9984 
summary(Choma5_clusts_shorter2_male$perc_time_in_all_clusts[!duplicated(Choma5_clusts_shorter2_male$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6629  0.8621  0.9165  0.8969  0.9323  0.9955 

# time outside clusters
summary(Choma5_clusts_shorter2_female$total_time[!duplicated(Choma5_clusts_shorter2_female$partid)] - Choma5_clusts_shorter2_female$part_trips_time_new[!duplicated(Choma5_clusts_shorter2_female$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1582  0.5242  1.3370  2.1695  3.8729  8.3344 
summary(Choma5_clusts_shorter2_male$total_time[!duplicated(Choma5_clusts_shorter2_male$partid)] - Choma5_clusts_shorter2_male$part_trips_time_new[!duplicated(Choma5_clusts_shorter2_male$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1395  1.6176  2.8666  3.3659  4.7557  9.8070 
#
### Time at locations ####
favstats(Choma5_clusts_shorter2_no_home_female$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2461111 1.155347 1.929583 5.299653 834.3897 35.29828 132.1338 322       0
favstats(Choma5_clusts_shorter2_no_home_male$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2605556 1.037222 1.937222 4.495903 885.6553 20.41479 96.68901 408       0

### PER PART ###
clust_time_part_median_female <- (favstats(Choma5_clusts_shorter2_no_home_female$clust_trips_time_new~Choma5_clusts_shorter2_no_home_female$partid)$median)*24
clust_time_part_median_male <- (favstats(Choma5_clusts_shorter2_no_home_male$clust_trips_time_new~Choma5_clusts_shorter2_no_home_male$partid)$median)*24
favstats(clust_time_part_median_female)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.8297222 1.344722 1.7825 2.410833 7.662778 2.373788 1.661407 29       0
favstats(clust_time_part_median_male)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.7801389 1.417639 1.923611 2.230069 3.972778 1.967509 0.7095319 31       0

clust_time_part_gender <- as.data.frame(cbind(c(clust_time_part_median_female, clust_time_part_median_male),
                                               c(rep("female",length(clust_time_part_median_female)),
                                                 rep("male",length(clust_time_part_median_male)))))
colnames(clust_time_part_gender) <- c("values","type")
clust_time_part_gender$values <- as.numeric(as.character(clust_time_part_gender$values))
clust_time_part_gender$type <- factor(clust_time_part_gender$type)
kruskal.test(clust_time_part_gender$values, clust_time_part_gender$type)

ggplot(data=clust_time_part_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,8)) +
  scale_x_continuous(breaks=c(seq(0,8,1)))
#
#
#
####### Biting time metrics #####
### number locations outside home ####
Choma5_clusts_shorter2_no_home_biting_places2_female <- Choma5_clusts_shorter2_no_home_biting_places2[which(Choma5_clusts_shorter2_no_home_biting_places2$male == 0),]
Choma5_clusts_shorter2_no_home_biting_places2_male <- Choma5_clusts_shorter2_no_home_biting_places2[which(Choma5_clusts_shorter2_no_home_biting_places2$male == 1),]

loc_counts_female <- favstats(Choma5_clusts_shorter2_no_home_biting_places2_female$clust~Choma5_clusts_shorter2_no_home_biting_places2_female$partid)$n
loc_counts_male <- favstats(Choma5_clusts_shorter2_no_home_biting_places2_male$clust~Choma5_clusts_shorter2_no_home_biting_places2_male$partid)$n
only_home_female <- nlevels(as.factor(Choma5_clusts_shorter2_female$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_female$partid))
only_home_male <- nlevels(as.factor(Choma5_clusts_shorter2_male$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_male$partid))
loc_counts_female2 <- c(loc_counts_female,rep(0, only_home_female))
loc_counts_male2 <- c(loc_counts_male,rep(0, only_home_male))

favstats(loc_counts_female2)
# min Q1 median Q3 max     mean       sd  n missing
#   0  0      1 2.75   6  1.7 1.765 30       0
favstats(loc_counts_male2)
# min Q1 median Q3 max     mean       sd  n missing
#   0 0.5      2  3   6 2.129 1.821 31       0
summary(as.factor(loc_counts_female2))
# 0 1 2 3 4 5 6 
# 9 8 5 4 1 1 2 
summary(as.factor(loc_counts_male2))
# 0 1 2 3 4 5 6 
# 8 6 3 7 3 3 1 
## all those with 0 have biting time at home location (not just nowhere)

locs_count2 <- as.data.frame(cbind(c(loc_counts_female2, loc_counts_male2),
                                   c(rep("female",length(loc_counts_female2)),
                                     rep("male",length(loc_counts_male2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-0.35,6.35)) +
  scale_x_continuous(breaks=c(seq(0,7,1)))
#
##
### time per location ####
favstats((Choma5_clusts_shorter2_no_home_biting_places2_female$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5636 3.052  16.22 137.5 291.3 63.44 87.91 51       0
favstats((Choma5_clusts_shorter2_no_home_biting_places2_male$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5142 0.9003  2.449 14.84 351.6 37.39 79.03 66       0

clusts_median_Choma_female <- (favstats(Choma5_clusts_shorter2_no_home_biting_places2_female$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_female$partid)$median)*24
clusts_median_Choma_male <- (favstats(Choma5_clusts_shorter2_no_home_biting_places2_male$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_male$partid)$median)*24
favstats(clusts_median_Choma_female)
# min       Q1   median  Q3      max     mean       sd  n missing
# 2.558 7.58  20.62 102.2 291.3 70.58 78.29 21       0
favstats(clusts_median_Choma_male)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.5928 1.784  2.733 10.93 244.5 36.31 74.9 23       0

clust_time_part_gender <- as.data.frame(cbind(c(clusts_median_Choma_female, clusts_median_Choma_male),
                                              c(rep("female",length(clusts_median_Choma_female)),
                                                rep("male",length(clusts_median_Choma_male)))))
colnames(clust_time_part_gender) <- c("values","type")
clust_time_part_gender$values <- as.numeric(as.character(clust_time_part_gender$values))
clust_time_part_gender$type <- factor(clust_time_part_gender$type)

ggplot(data=clust_time_part_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=2.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,300)) +
  scale_x_continuous(breaks=c(seq(0,300,50)))
#
### number trips ####
Choma5_clust_trips_shorter2_no_home_biting_places3_female <- Choma5_clust_trips_shorter2_no_home_biting_places3[which(Choma5_clust_trips_shorter2_no_home_biting_places3$male == 0),]
Choma5_clust_trips_shorter2_no_home_biting_places3_male <- Choma5_clust_trips_shorter2_no_home_biting_places3[which(Choma5_clust_trips_shorter2_no_home_biting_places3$male == 1),]

favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#    1  1      2 3.5  23 4.294 5.608 51       0
favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  3  27 3.969 6.359 65       0

summary(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2))
#  1  2  3  4  6  7  8 10 12 15 22 23 
# 18 16  4  1  1  1  2  3  1  1  2  1 
summary(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2))
#  1  2  3  4  5  8 12 13 14 17 21 22 27 
# 43  3  4  2  3  1  1  2  1  1  1  1  2 

### PER PART ###
trip_count_median_female <- favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_female$partid)$median
favstats(trip_count_median_female)
# min Q1 median Q3 max  mean  sd  n missing
#   1  2      2 4.5  12 3.31 2.62 21       0
trip_count_median_male <- favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_male$partid)$median
favstats(trip_count_median_male)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1 2.5  27 3.957 6.807 23       0

trips_part2_gender <- as.data.frame(cbind(c(trip_count_median_female, trip_count_median_male),
                                          c(rep("female",length(trip_count_median_female)),
                                            rep("male",length(trip_count_median_male)))))
colnames(trips_part2_gender) <- c("values","type")
trips_part2_gender$values <- as.numeric(as.character(trips_part2_gender$values))
trips_part2_gender$type <- factor(trips_part2_gender$type)

ggplot(data=trips_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,26.25)) + 
  scale_x_continuous(breaks=c(seq(0,30,2)))


## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,30))+ 
  scale_x_continuous(breaks=c(seq(0,30,2)))
#
##### percent biting time spent at home vs elsewhere #####
Choma5_clusts_shorter2_biting_places3b_female <- Choma5_clusts_shorter2_biting_places3b[which(Choma5_clusts_shorter2_biting_places3b$male == 0),]
Choma5_clusts_shorter2_biting_places3b_male <- Choma5_clusts_shorter2_biting_places3b[which(Choma5_clusts_shorter2_biting_places3b$male == 1),]

favstats(Choma5_clusts_shorter2_biting_places3b_female$percent_home_biting_time)
favstats(Choma5_clusts_shorter2_biting_places3b_female$percent_home_biting_time[which(Choma5_clusts_shorter2_biting_places3b_female$percent_home_biting_time != 0)])
favstats(Choma5_clusts_shorter2_biting_places3b_male$percent_home_biting_time)
favstats(Choma5_clusts_shorter2_biting_places3b_male$percent_home_biting_time[which(Choma5_clusts_shorter2_biting_places3b_male$percent_home_biting_time != 0)])


perc_time_home2_gender <- as.data.frame(cbind(c(Choma5_clusts_shorter2_biting_places3b_female$percent_home_biting_time, Choma5_clusts_shorter2_biting_places3b_male$percent_home_biting_time),
                                              c(rep("female",length(Choma5_clusts_shorter2_biting_places3b_female$percent_home_biting_time)),
                                                rep("male",length(Choma5_clusts_shorter2_biting_places3b_male$percent_home_biting_time)))))
colnames(perc_time_home2_gender) <- c("values","type")
perc_time_home2_gender$values <- as.numeric(as.character(perc_time_home2_gender$values))
perc_time_home2_gender$type <- factor(perc_time_home2_gender$type)

ggplot(data=perc_time_home2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,100))

ggplot(data=Choma5_clusts_shorter2_biting_places3b, aes(x=(percent_home_biting_time))) +
  geom_histogram(binwidth = 10, fill = "black") +
  coord_cartesian(xlim=c(0,100)) +
  scale_x_continuous(breaks=c(seq(0,100,10)))
##
###############
###############
######### MUTASA #########
###### General Gender Analysis ####
Mutasa_part_info <- Mutasa5_clusts_shorter2[!duplicated(Mutasa5_clusts_shorter2$partid),c(1:3)]
summary(as.factor(Mutasa_part_info$male)) 
#    0    1
#  94    88
Mutasa_part_info$partid2 <- gsub("b","", Mutasa_part_info$partid)
Mutasa_part_info2 <- Mutasa_part_info[!duplicated(Mutasa_part_info$partid2),]
summary(as.factor(Mutasa_part_info2$male))
###### Split by Gender ######
### number of locations (without home) ####
Mutasa5_clusts_shorter2_no_home_female <- Mutasa5_clusts_shorter2_no_home[which(Mutasa5_clusts_shorter2_no_home$male == 0),]
Mutasa5_clusts_shorter2_no_home_male <- Mutasa5_clusts_shorter2_no_home[which(Mutasa5_clusts_shorter2_no_home$male == 1),]
loc_counts_female <- favstats(Mutasa5_clusts_shorter2_no_home_female$clust~Mutasa5_clusts_shorter2_no_home_female$partid)$n
loc_counts_male <- favstats(Mutasa5_clusts_shorter2_no_home_male$clust~Mutasa5_clusts_shorter2_no_home_male$partid)$n
only_home_female <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_female$partid))
only_home_male <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_male$partid))
loc_counts_female2 <- c(loc_counts_female,rep(0, only_home_female))
loc_counts_male2 <- c(loc_counts_male,rep(0, only_home_male))
favstats(loc_counts_female2)
# min   Q1 median     Q3  max     mean       sd  n missing
# 0  3      6 9.75  24 6.87234 4.980033 94       0
favstats(loc_counts_male2)
# min Q1 median Q3 max     mean       sd  n missing
#    0  4      8 13.75  25 9.104651 6.264802 86       0

locs_count2 <- as.data.frame(cbind(c(loc_counts_female2, loc_counts_male2),
                                   c(rep("female",length(loc_counts_female2)),
                                     rep("male",length(loc_counts_male2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

kruskal.test(locs_count2$values, locs_count2$type) # p = 0.02

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-2,27)) +
  scale_x_continuous(breaks=c(seq(0,25,5)))


my_comp <- list(c("female","male"))
ggplot(data=locs_count2, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("red","blue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 30) 

#
### distance of locations (without home) ####
### TOTAL ###
favstats(Mutasa5_clusts_shorter2_no_home_female$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#  2.745604e-05 0.4424186 1.23664 3.871323 228.4774 13.82587 44.09112 646       0
favstats(Mutasa5_clusts_shorter2_no_home_male$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 4.190736e-05 0.5494893 1.409165 4.736531 481.4493 19.72159 62.22123 783       0

### PER PART ###
medians_km_female <- (favstats(Mutasa5_clusts_shorter2_no_home_female$hhdist_m_haversine~Mutasa5_clusts_shorter2_no_home_female$partid)$median)/1000
medians_km_male <- (favstats(Mutasa5_clusts_shorter2_no_home_male$hhdist_m_haversine~Mutasa5_clusts_shorter2_no_home_male$partid)$median)/1000
favstats(medians_km_female) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
# 0.08078032 0.5716804 1.234345 2.716702 217.5585 9.874464 34.74559 91       0
favstats(medians_km_male) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
# 0.1108818 0.7020393 1.072867 2.492948 299.3132 11.37001 40.76159 83       0

## REMOVE OUTLIER ##
medians_km_female_short <- medians_km_female[-which(medians_km_female > 100)]
favstats(medians_km_female_short)
# min        Q1    median       Q3      max     mean       sd  n missing
# 0.08078032 0.5670885 1.190954 2.238006 36.2448 2.939855 5.648665 87       0

medians_km_male_short <- medians_km_male[-which(medians_km_male > 100)]
favstats(medians_km_male_short)
# min        Q1    median       Q3      max     mean       sd  n missing
# 0.1108818 0.7013227 1.06786 2.276097 57.42812 5.365199 11.7424 81       0

dists_km2_gender <- as.data.frame(cbind(c(medians_km_female_short, medians_km_male_short),
                                        c(rep("female",length(medians_km_female_short)),
                                          rep("male",length(medians_km_male_short)))))
colnames(dists_km2_gender) <- c("values","type")
dists_km2_gender$values <- as.numeric(as.character(dists_km2_gender$values))
dists_km2_gender$type <- factor(dists_km2_gender$type)
kruskal.test(dists_km2_gender$values, dists_km2_gender$type) 

ggplot(data=dists_km2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=10) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-2,65))
#
### number of trips (without home) ####
### TOTAL ###
Mutasa5_clust_trips_shorter2_no_home3_female <- Mutasa5_clust_trips_shorter2_no_home3[which(Mutasa5_clust_trips_shorter2_no_home3$male == 0),]
Mutasa5_clust_trips_shorter2_no_home3_male <- Mutasa5_clust_trips_shorter2_no_home3[which(Mutasa5_clust_trips_shorter2_no_home3$male == 1),]

favstats(Mutasa5_clust_trips_shorter2_no_home3_female$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      1  2  32 2.140867 2.412558 646       0
favstats(Mutasa5_clust_trips_shorter2_no_home3_male$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      2  3  43 2.913155 4.246571 783       0

### PER PART ###
trips_part_median_female <- favstats(Mutasa5_clust_trips_shorter2_no_home3_female$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_female$partid)$median
trips_part_median_male <- favstats(Mutasa5_clust_trips_shorter2_no_home3_male$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_male$partid)$median
favstats(trips_part_median_female)
#   min  Q1   median   Q3      max       mean       sd  n missing
#   1  1      1  2  13 1.686813 1.514097 91       0
favstats(trips_part_median_male)
#   min Q1 median  Q3   max      mean       sd  n missing
#   1  1      2  2  13 2.120482 1.812205 83       0

trips_part2_gender <- as.data.frame(cbind(c(trips_part_median_female, trips_part_median_male),
                                          c(rep("female",length(trips_part_median_female)),
                                            rep("male",length(trips_part_median_male)))))
colnames(trips_part2_gender) <- c("values","type")
trips_part2_gender$values <- as.numeric(as.character(trips_part2_gender$values))
trips_part2_gender$type <- factor(trips_part2_gender$type)

kruskal.test(trips_part2_gender$values, trips_part2_gender$type) # p = 0.002

ggplot(data=trips_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,15)) +
  scale_x_continuous(breaks=c(0:15))
#
### time per trip ####
### TOTAL ###
Mutasa5_clust_trips_shorter2_no_home_female <- Mutasa5_clust_trips_shorter2_no_home[which(Mutasa5_clust_trips_shorter2_no_home$male == 0),]
Mutasa5_clust_trips_shorter2_no_home_male <- Mutasa5_clust_trips_shorter2_no_home[which(Mutasa5_clust_trips_shorter2_no_home$male == 1),]
favstats(Mutasa5_clust_trips_shorter2_no_home_female$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
#0.4413889 1.639167 3.067222 5.76 856.8375 8.063954 33.01965 1383       0
favstats(Mutasa5_clust_trips_shorter2_no_home_male$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.48 1.270278 2.367778 5.085833 1218.114 6.260339 35.98856 2281       0

### PER PART ###
trip_time_part_median_female <- favstats(Mutasa5_clust_trips_shorter2_no_home_female$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_female$partid)$median
trip_time_part_median_male <- favstats(Mutasa5_clust_trips_shorter2_no_home_male$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_male$partid)$median
favstats(trip_time_part_median_female)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.9602778 2.254722 3.049722 3.965903 33.96292 3.865583 3.823232 91       0
favstats(trip_time_part_median_male)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.066667 1.984931 2.498333 3.375278 1218.114 17.83338 133.3692 83       0

## REMOVE OUTLIER ##
trip_time_part_median_female_short <- trip_time_part_median_female[-which(trip_time_part_median_female > 30)]
favstats(trip_time_part_median_female_short)
# min       Q1   median       Q3      max     mean       sd   n missing
# 0.9602778 2.24625 3.018889 3.934306 12.22667 3.531168 2.119001 90       0

trip_time_part_median_male_short <- trip_time_part_median_male[-which(trip_time_part_median_male > 600)]
favstats(trip_time_part_median_male_short)
# min       Q1   median       Q3      max     mean       sd   n missing
# 1.066667 1.984271 2.487361 3.330347 9.897222 3.195815 1.986377 82       0

trip_times_part2_gender <- as.data.frame(cbind(c(trip_time_part_median_female_short, trip_time_part_median_male_short),
                                               c(rep("female",length(trip_time_part_median_female_short)),
                                                 rep("male",length(trip_time_part_median_male_short)))))
colnames(trip_times_part2_gender) <- c("values","type")
trip_times_part2_gender$values <- as.numeric(as.character(trip_times_part2_gender$values))
trip_times_part2_gender$type <- factor(trip_times_part2_gender$type)
kruskal.test(trip_times_part2_gender$values, trip_times_part2_gender$type) 

ggplot(data=trip_times_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-1,13)) +
  scale_x_continuous(breaks=c(seq(0,15,2)))

#
### Home Loc ####
Mutasa5_clusts_shorter2_home_female <- Mutasa5_clusts_shorter2_home[which(Mutasa5_clusts_shorter2_home$male == 0),] 
Mutasa5_clusts_shorter2_home_male <- Mutasa5_clusts_shorter2_home[which(Mutasa5_clusts_shorter2_home$male == 1),] 

no_home_female <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_home_female$partid))
## no females with no home
no_home_male <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_home_male$partid))
## no males with no homes 

### number of 'trips' to home loc ####
### TOTAL ###
Mutasa5_clust_trips_shorter2_home3_female <- Mutasa5_clust_trips_shorter2_home3[which(Mutasa5_clust_trips_shorter2_home3$male == 0),]
Mutasa5_clust_trips_shorter2_home3_male <- Mutasa5_clust_trips_shorter2_home3[which(Mutasa5_clust_trips_shorter2_home3$male == 1),]

favstats(Mutasa5_clust_trips_shorter2_home3_female$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#    1  5      8 15  31 10.51613 7.517534 93       0
favstats(Mutasa5_clust_trips_shorter2_home3_male$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#     1  6     12 23  42 14.62353 11.0658 85       0

### total time at home loc ####
Mutasa5_clust_trips_shorter2_home4_female <- Mutasa5_clust_trips_shorter2_home4[which(Mutasa5_clust_trips_shorter2_home4$male == 0),]
Mutasa5_clust_trips_shorter2_home4_male <- Mutasa5_clust_trips_shorter2_home4[which(Mutasa5_clust_trips_shorter2_home4$male == 1),]

favstats(Mutasa5_clust_trips_shorter2_home4_female$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
# 0.920081 19.14721 25.64502 32.52499 59.14673 26.43902 11.03208 93       0
favstats(Mutasa5_clust_trips_shorter2_home4_male$clust_trips_time_new)
#         min       Q1   median       Q3      max     mean       sd  n missing
# 0.02479167 14.443 21.37396 27.03878 65.96275 22.57999 11.98395 85       0
favstats(Mutasa5_clust_trips_shorter2_home4_male$clust_trips_time_new[which(Mutasa5_clust_trips_shorter2_home4_male$partid != "30240")])

### percent time at home ####
favstats(Mutasa5_clusts_shorter2_home_female$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
#   14.3927 75.75095 91.9012 97.67469 100 85.52081 16.72762 94       0
favstats(Mutasa5_clusts_shorter2_home_male$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#   16.36831 69.6522 82.56211 93.53527 100 77.71848 19.99931 86       0

perc_time_home2_gender <- as.data.frame(cbind(c(Mutasa5_clusts_shorter2_home_female$perc_time_home, Mutasa5_clusts_shorter2_home_male$perc_time_home),
                                              c(rep("female",length(Mutasa5_clusts_shorter2_home_female$perc_time_home)),
                                                rep("male",length(Mutasa5_clusts_shorter2_home_male$perc_time_home)))))
colnames(perc_time_home2_gender) <- c("values","type")
perc_time_home2_gender$values <- as.numeric(as.character(perc_time_home2_gender$values))
perc_time_home2_gender$type <- factor(perc_time_home2_gender$type)

kruskal.test(perc_time_home2_gender$values, perc_time_home2_gender$type) # p = 0.002

ggplot(data=perc_time_home2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,100))
##
my_comp <- list(c("female","male"))
ggplot(data=perc_time_home2_gender, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("red","blue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 115) 
#
##
### Average Direction of Travel #####
Mutasa5_clusts_shorter2_no_home$male <- as.factor(Mutasa5_clusts_shorter2_no_home$male)
temp <- Mutasa5_clusts_shorter2_no_home[!duplicated(Mutasa5_clusts_shorter2_no_home$partid),]
temp2 <- temp[which(!(is.na(temp$male))),]

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=male), position=position_dodge2(preserve = "single"), width=3) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Male", values=c("deeppink","blue"))
#

favstats(temp2$avg_dir_from_home_bearing~temp2$male)
# temp2$male       min       Q1   median       Q3      max     mean        sd  n missing
# 1          0 16.387774 137.9740 221.8256 290.5334 359.1797 207.7788  99.14697 91       0
# 2          1  2.636555 117.7448 221.1990 276.4839 359.7089 196.4702 102.04706 83       0
favstats(temp2$avg_dir_from_home_mag~temp2$male)
# temp2$male        min        Q1   median       Q3      max     mean       sd  n missing
# 1          0 0.03380367 0.6042695 1.674874 5.700955 210.0080 12.08885 33.33847 91       0
# 2          1 0.09139251 0.6000253 1.244584 5.019252 181.9573 15.05311 36.58184 83       0

temp3 <- temp2[which(temp2$avg_dir_from_home_mag < 70),] ## 164 of 176 parts

favstats(temp3$avg_dir_from_home_bearing~temp3$male)
# temp3$male       min       Q1   median       Q3      max     mean        sd  n missing
# 1          0 16.387774 126.2034 219.9242 289.0636 359.1797 204.1143  99.88831 87       0
# 2          1  2.636555 112.8585 211.7342 272.4138 359.7089 192.4088 102.92370 75       0
favstats(temp3$avg_dir_from_home_mag~temp3$male)
# temp3$male        min        Q1   median       Q3      max     mean        sd  n missing
# 1          0 0.03380367 0.5612283 1.437137 4.436135 65.41056 5.854193 11.909021 87       0
# 2          1 0.09139251 0.5655562 1.182748 3.112735 68.50646 4.357640  9.733805 75       0
kruskal.test(temp3$avg_dir_from_home_bearing, temp3$male) 

ggplot(data=temp3, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=male), position=position_dodge2(preserve = "single"), width=3) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Male", values=c("red","blue"))





temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)] <- -(360-temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)])
ggplot(data=temp2, aes(x=avg_dir_from_home_values, y=..count.., fill=male)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~male, labeller = labeller(male = c("0"="Female", "1"="Male")))+
  scale_fill_manual(values=c("red","blue")) +
  scale_y_continuous(breaks=c(seq(0,14,2)))

#
#

### time spent at clusters vs elsewhere (travel, less important places, etc) ####
Mutasa5_clusts_shorter2_female <- Mutasa5_clusts_shorter2[which(Mutasa5_clusts_shorter2$male == 0),]
Mutasa5_clusts_shorter2_male <- Mutasa5_clusts_shorter2[which(Mutasa5_clusts_shorter2$male == 1),]
Mutasa5_clusts_shorter2_female$perc_time_in_all_clusts <- Mutasa5_clusts_shorter2_female$part_trips_time_new/Mutasa5_clusts_shorter2_female$total_time
Mutasa5_clusts_shorter2_male$perc_time_in_all_clusts <- Mutasa5_clusts_shorter2_male$part_trips_time_new/Mutasa5_clusts_shorter2_male$total_time

summary(Mutasa5_clusts_shorter2_female$perc_time_in_all_clusts[!duplicated(Mutasa5_clusts_shorter2_female$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7556  0.8992  0.9560  0.9374  0.9872  0.9984 
summary(Mutasa5_clusts_shorter2_male$perc_time_in_all_clusts[!duplicated(Mutasa5_clusts_shorter2_male$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6629  0.8621  0.9165  0.8969  0.9323  0.9955 

# time outside clusters
summary(Mutasa5_clusts_shorter2_female$total_time[!duplicated(Mutasa5_clusts_shorter2_female$partid)] - Mutasa5_clusts_shorter2_female$part_trips_time_new[!duplicated(Mutasa5_clusts_shorter2_female$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1582  0.5242  1.3370  2.1695  3.8729  8.3344 
summary(Mutasa5_clusts_shorter2_male$total_time[!duplicated(Mutasa5_clusts_shorter2_male$partid)] - Mutasa5_clusts_shorter2_male$part_trips_time_new[!duplicated(Mutasa5_clusts_shorter2_male$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1395  1.6176  2.8666  3.3659  4.7557  9.8070 

## part 30320b --> ~1/3 of points considered 'outliers' --> seem to be around home, but spread out enough to not count
#
### Time at locations ####
favstats(Mutasa5_clusts_shorter2_no_home_female$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.4802778 2.342917 4.233472 8.631597 856.8375 17.26385 64.75601 646       0
favstats(Mutasa5_clusts_shorter2_no_home_male$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.5072222 2.227222 4.079722 9.601944 1537.291 18.23733 81.05983 783       0

### PER PART ###
clust_time_part_median_female <- (favstats(Mutasa5_clusts_shorter2_no_home_female$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_female$partid)$median)*24
clust_time_part_median_male <- (favstats(Mutasa5_clusts_shorter2_no_home_male$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_male$partid)$median)*24
favstats(clust_time_part_median_female)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 1.656389 3.208264 3.960278 5.792986 98.0725 6.81656 12.55591 91       0
favstats(clust_time_part_median_male)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.433611 3.354722 4.354861 6.473472 1218.114 24.87472 135.0354 83       0
clust_time_part_median_male2 <- clust_time_part_median_male[which(clust_time_part_median_male < 150)]
favstats(clust_time_part_median_male2)
# min       Q1   median       Q3      max     mean       sd  n missing
# 1.433611 3.349444 4.343333 6.416389 110.4386 8.038273 15.56373 81       0

clust_time_part_gender <- as.data.frame(cbind(c(clust_time_part_median_female, clust_time_part_median_male2),
                                              c(rep("female",length(clust_time_part_median_female)),
                                                rep("male",length(clust_time_part_median_male2)))))
colnames(clust_time_part_gender) <- c("values","type")
clust_time_part_gender$values <- as.numeric(as.character(clust_time_part_gender$values))
clust_time_part_gender$type <- factor(clust_time_part_gender$type)
kruskal.test(clust_time_part_gender$values, clust_time_part_gender$type) 

ggplot(data=clust_time_part_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-1,115)) +
  scale_x_continuous(breaks=c(seq(0,120,10)))
#
#
#
####### Biting time metrics #####
### number locations outside home ####
Mutasa5_clusts_shorter2_no_home_biting_places2_female <- Mutasa5_clusts_shorter2_no_home_biting_places2[which(Mutasa5_clusts_shorter2_no_home_biting_places2$male == 0),]
Mutasa5_clusts_shorter2_no_home_biting_places2_male <- Mutasa5_clusts_shorter2_no_home_biting_places2[which(Mutasa5_clusts_shorter2_no_home_biting_places2$male == 1),]

loc_counts_female <- favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_female$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_female$partid)$n
loc_counts_male <- favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_male$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_male$partid)$n
only_home_female <- nlevels(as.factor(Mutasa5_clusts_shorter2_female$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_female$partid))
only_home_male <- nlevels(as.factor(Mutasa5_clusts_shorter2_male$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_male$partid))
loc_counts_female2 <- c(loc_counts_female,rep(0, only_home_female))
loc_counts_male2 <- c(loc_counts_male,rep(0, only_home_male))

favstats(loc_counts_female2)
# min Q1 median Q3 max     mean       sd  n missing
#   0  0    0.5  2   5 0.9894 1.231 94       0
favstats(loc_counts_male2)
# min Q1 median Q3 max     mean       sd  n missing
#  0  0      1 2.75   8 1.849 2.134 86       0
summary(as.factor(loc_counts_female2))
#  0  1  2   3 4  5 
# 47 18 18  6  4  1 
summary(as.factor(loc_counts_male2))
#  0  1  2  3  4  6  7  8 
# 26 24 14  8  5  4  1  4 
## all those with 0 have biting time at home location (not just nowhere)

locs_count2 <- as.data.frame(cbind(c(loc_counts_female2, loc_counts_male2),
                                   c(rep("female",length(loc_counts_female2)),
                                     rep("male",length(loc_counts_male2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-0.35,8.35)) +
  scale_x_continuous(breaks=c(seq(0,8,1)))
#
##
### time per location ####
favstats((Mutasa5_clusts_shorter2_no_home_biting_places2_female$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
#  0.5072 4.096  8.752 17.88 264.9 24.65 44.6 93       0
favstats((Mutasa5_clusts_shorter2_no_home_biting_places2_male$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5222 1.63  5.528 12.43 504.1 19.39 49.81 159       0

clusts_median_Mutasa_female <- (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_female$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_female$partid)$median)*24
clusts_median_Mutasa_male <- (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_male$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_male$partid)$median)*24
favstats(clusts_median_Mutasa_female)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.5946 5.682  8.936 30.67 234.9 25.37 40.77 47       0
favstats(clusts_median_Mutasa_male)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.5831 2.624  7.272 10.71 504.1 24.25 69.9 60       0

clust_time_part_gender <- as.data.frame(cbind(c(clusts_median_Mutasa_female, clusts_median_Mutasa_male),
                                              c(rep("female",length(clusts_median_Mutasa_female)),
                                                rep("male",length(clusts_median_Mutasa_male)))))
colnames(clust_time_part_gender) <- c("values","type")
clust_time_part_gender$values <- as.numeric(as.character(clust_time_part_gender$values))
clust_time_part_gender$type <- factor(clust_time_part_gender$type)

ggplot(data=clust_time_part_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=2.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,500)) +
  scale_x_continuous(breaks=c(seq(0,500,50)))
#

clust_time_part_gender2 <- clust_time_part_gender[which(clust_time_part_gender$values < 500),]
ggplot(data=clust_time_part_gender2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=2.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,250)) +
  scale_x_continuous(breaks=c(seq(0,300,50)))
#


### number trips ####
Mutasa5_clust_trips_shorter2_no_home_biting_places3_female <- Mutasa5_clust_trips_shorter2_no_home_biting_places3[which(Mutasa5_clust_trips_shorter2_no_home_biting_places3$male == 0),]
Mutasa5_clust_trips_shorter2_no_home_biting_places3_male <- Mutasa5_clust_trips_shorter2_no_home_biting_places3[which(Mutasa5_clust_trips_shorter2_no_home_biting_places3$male == 1),]

favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2  14 1.806 1.969 93       0
favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#  1  1      1  2  22 2.421 3.65 159       0

summary(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2))
#  1  2  3  4  5  6  7  8  9 14 
# 66 12  7  2  1  1  1  1  1  1 
summary(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2))
#   1   2   3   4   6   7   8   9  10  11  13  14  22 
# 108  24   6   4   2   1   3   3   1   2   1   1   3 

### PER PART ###
trip_count_median_female <- favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_female$partid)$median
favstats(trip_count_median_female)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1 1.5   6 1.532 1.115 47       0
trip_count_median_male <- favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_male$partid)$median
favstats(trip_count_median_male)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2  13 2.025 2.438 60       0

trips_part2_gender <- as.data.frame(cbind(c(trip_count_median_female, trip_count_median_male),
                                          c(rep("female",length(trip_count_median_female)),
                                            rep("male",length(trip_count_median_male)))))
colnames(trips_part2_gender) <- c("values","type")
trips_part2_gender$values <- as.numeric(as.character(trips_part2_gender$values))
trips_part2_gender$type <- factor(trips_part2_gender$type)

ggplot(data=trips_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,13)) + 
  scale_x_continuous(breaks=c(seq(0,30,2)))


## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,30))+ 
  scale_x_continuous(breaks=c(seq(0,30,2)))
#
##### percent biting time spent at home vs elsewhere #####
Mutasa5_clusts_shorter2_biting_places3b_female <- Mutasa5_clusts_shorter2_biting_places3b[which(Mutasa5_clusts_shorter2_biting_places3b$male == 0),]
Mutasa5_clusts_shorter2_biting_places3b_male <- Mutasa5_clusts_shorter2_biting_places3b[which(Mutasa5_clusts_shorter2_biting_places3b$male == 1),]

favstats(Mutasa5_clusts_shorter2_biting_places3b_female$percent_home_biting_time)
favstats(Mutasa5_clusts_shorter2_biting_places3b_female$percent_home_biting_time[which(Mutasa5_clusts_shorter2_biting_places3b_female$percent_home_biting_time != 0)])
favstats(Mutasa5_clusts_shorter2_biting_places3b_male$percent_home_biting_time)
favstats(Mutasa5_clusts_shorter2_biting_places3b_male$percent_home_biting_time[which(Mutasa5_clusts_shorter2_biting_places3b_male$percent_home_biting_time != 0)])


perc_time_home2_gender <- as.data.frame(cbind(c(Mutasa5_clusts_shorter2_biting_places3b_female$percent_home_biting_time, Mutasa5_clusts_shorter2_biting_places3b_male$percent_home_biting_time),
                                              c(rep("female",length(Mutasa5_clusts_shorter2_biting_places3b_female$percent_home_biting_time)),
                                                rep("male",length(Mutasa5_clusts_shorter2_biting_places3b_male$percent_home_biting_time)))))
colnames(perc_time_home2_gender) <- c("values","type")
perc_time_home2_gender$values <- as.numeric(as.character(perc_time_home2_gender$values))
perc_time_home2_gender$type <- factor(perc_time_home2_gender$type)

ggplot(data=perc_time_home2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,100))
##
###############
###############
##### NCHELENGE #####
###### General Gender Analysis ####
Nchelenge_part_info <- Nchelenge5_clusts_shorter2[!duplicated(Nchelenge5_clusts_shorter2$partid),c(1:3)]
summary(as.factor(Nchelenge_part_info$male)) 
#   0   1
#  45  30
###### Split by Gender ######
### number of locations (without home) ####
Nchelenge5_clusts_shorter2_no_home_female <- Nchelenge5_clusts_shorter2_no_home[which(Nchelenge5_clusts_shorter2_no_home$male == 0),]
Nchelenge5_clusts_shorter2_no_home_male <- Nchelenge5_clusts_shorter2_no_home[which(Nchelenge5_clusts_shorter2_no_home$male == 1),]
loc_counts_female <- favstats(Nchelenge5_clusts_shorter2_no_home_female$clust~Nchelenge5_clusts_shorter2_no_home_female$partid)$n
loc_counts_male <- favstats(Nchelenge5_clusts_shorter2_no_home_male$clust~Nchelenge5_clusts_shorter2_no_home_male$partid)$n
only_home_female <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_female$partid))
only_home_male <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_male$partid))
loc_counts_female2 <- c(loc_counts_female,rep(0, only_home_female))
loc_counts_male2 <- c(loc_counts_male,rep(0, only_home_male))
favstats(loc_counts_female2)
# min   Q1 median     Q3  max     mean       sd  n missing
#    0  3      8 12  19 8.177778 5.453476 45       0
favstats(loc_counts_male2)
# min Q1 median Q3 max     mean       sd  n missing
#   0 5.25    8.5 12  23 9.166667 6.069047 30       0

locs_count2 <- as.data.frame(cbind(c(loc_counts_female2, loc_counts_male2),
                                   c(rep("female",length(loc_counts_female2)),
                                     rep("male",length(loc_counts_male2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

kruskal.test(locs_count2$values, locs_count2$type) # p = 0.53

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-2,27)) +
  scale_x_continuous(breaks=c(seq(0,25,5)))

#
### distance of locations (without home) ####
### TOTAL ###
favstats(Nchelenge5_clusts_shorter2_no_home_female$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#  0 0.3613393 1.030946 3.325195 187.2311 4.927711 15.93934 368       0
favstats(Nchelenge5_clusts_shorter2_no_home_male$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 0.02460256 0.3722302 1.497811 7.836124 146.5504 9.419034 24.45845 275       0

### PER PART ###
medians_km_female <- (favstats(Nchelenge5_clusts_shorter2_no_home_female$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_female$partid)$median)
medians_km_male <- (favstats(Nchelenge5_clusts_shorter2_no_home_male$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_male$partid)$median)
favstats(medians_km_female) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
# 0.1811803 0.5990128 1.185105 2.106645 186.7748 6.961053 28.78531 42       0
favstats(medians_km_male) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.1012507 0.5408554 1.12657 3.353959 85.42321 6.159835 16.20602 29       0

## REMOVE OUTLIER ##
medians_km_female_short <- medians_km_female[-which(medians_km_female > 70)]
favstats(medians_km_female_short)
# min        Q1    median       Q3      max     mean       sd  n missing
# 0.1811803 0.5938694 1.166392 2.046623 25.71711 2.575353 4.611859 41       0

medians_km_male_short <- medians_km_male[-which(medians_km_male > 70)]
favstats(medians_km_male_short)
# min        Q1    median       Q3      max     mean       sd  n missing
# 0.1012507 0.539771 1.045539 3.214469 24.49658 3.329 5.599993 28       0

dists_km2_gender <- as.data.frame(cbind(c(medians_km_female_short, medians_km_male_short),
                                        c(rep("female",length(medians_km_female_short)),
                                          rep("male",length(medians_km_male_short)))))
colnames(dists_km2_gender) <- c("values","type")
dists_km2_gender$values <- as.numeric(as.character(dists_km2_gender$values))
dists_km2_gender$type <- factor(dists_km2_gender$type)

kruskal.test(dists_km2_gender$values, dists_km2_gender$type) # p = 0.97

ggplot(data=dists_km2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,30))



dists_km2_gender <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_no_home_female$hhdist_km_haversine, Nchelenge5_clusts_shorter2_no_home_male$hhdist_km_haversine),
                                        c(rep("female",length(Nchelenge5_clusts_shorter2_no_home_female$hhdist_km_haversine)),
                                          rep("male",length(Nchelenge5_clusts_shorter2_no_home_male$hhdist_km_haversine)))))
colnames(dists_km2_gender) <- c("values","type")
dists_km2_gender$values <- as.numeric(as.character(dists_km2_gender$values))
dists_km2_gender$type <- factor(dists_km2_gender$type)

kruskal.test(dists_km2_gender$values, dists_km2_gender$type) # p = 0.016

ggplot(data=dists_km2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,30))
#
### number of trips (without home) ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2_no_home3_female <- Nchelenge5_clust_trips_shorter2_no_home3[which(Nchelenge5_clust_trips_shorter2_no_home3$male == 0),]
Nchelenge5_clust_trips_shorter2_no_home3_male <- Nchelenge5_clust_trips_shorter2_no_home3[which(Nchelenge5_clust_trips_shorter2_no_home3$male == 1),]

favstats(Nchelenge5_clust_trips_shorter2_no_home3_female$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      2  2  23 2.3125 2.454525 368       0
favstats(Nchelenge5_clust_trips_shorter2_no_home3_male$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      2  3  20 2.786232 3.005684 276       0

### PER PART ###
trips_part_median_female <- favstats(Nchelenge5_clust_trips_shorter2_no_home3_female$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_female$partid)$median
trips_part_median_male <- favstats(Nchelenge5_clust_trips_shorter2_no_home3_male$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_male$partid)$median
favstats(trips_part_median_female)
#    min Q1 median  Q3   max      mean       sd  n missing
#   1  1   1.75  2   3 1.619048 0.5157226 42       0
favstats(trips_part_median_male)
#   min Q1 median  Q3   max      mean       sd  n missing
#     1  1      2  2   4 1.672414 0.7473269 29       0

trips_part2_gender <- as.data.frame(cbind(c(trips_part_median_female, trips_part_median_male),
                                          c(rep("female",length(trips_part_median_female)),
                                            rep("male",length(trips_part_median_male)))))
colnames(trips_part2_gender) <- c("values","type")
trips_part2_gender$values <- as.numeric(as.character(trips_part2_gender$values))
trips_part2_gender$type <- factor(trips_part2_gender$type)

kruskal.test(trips_part2_gender$values, trips_part2_gender$type) # p = 1

ggplot(data=trips_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks=c(0:10))
#
### time per trip ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2_no_home_female <- Nchelenge5_clust_trips_shorter2_no_home[which(Nchelenge5_clust_trips_shorter2_no_home$male == 0),]
Nchelenge5_clust_trips_shorter2_no_home_male <- Nchelenge5_clust_trips_shorter2_no_home[which(Nchelenge5_clust_trips_shorter2_no_home$male == 1),]
favstats(Nchelenge5_clust_trips_shorter2_no_home_female$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2483333 0.59625 1.244167 3.110833 394.7653 7.192904 30.20011 851       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_male$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2416667 0.5456944 1.126667 2.303333 244.0533 3.982926 14.17849 759       0

### PER PART ###
trip_time_part_median_female <- favstats(Nchelenge5_clust_trips_shorter2_no_home_female$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_female$partid)$median
trip_time_part_median_male <- favstats(Nchelenge5_clust_trips_shorter2_no_home_male$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_male$partid)$median
favstats(trip_time_part_median_female)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4844444 0.9635417 1.436806 1.997292 6.148333 1.611124 0.9514975 42       0
favstats(trip_time_part_median_male)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.4955556 0.7780556 1.126667 1.68625 6.663194 1.533329 1.321883 29       0

## REMOVE OUTLIERS ##
trip_time_part_median_female_short <- trip_time_part_median_female ## no outliers to remove

trip_times_part2_gender <- as.data.frame(cbind(c(trip_time_part_median_female_short, trip_time_part_median_male),
                                               c(rep("female",length(trip_time_part_median_female_short)),
                                                 rep("male",length(trip_time_part_median_male)))))
colnames(trip_times_part2_gender) <- c("values","type")
trip_times_part2_gender$values <- as.numeric(as.character(trip_times_part2_gender$values))
trip_times_part2_gender$type <- factor(trip_times_part2_gender$type)
kruskal.test(trip_times_part2_gender$values, trip_times_part2_gender$type) 

ggplot(data=trip_times_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 0.5, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,7)) +
  scale_x_continuous(breaks=c(0:7))
#
### Home Loc ####
Nchelenge5_clusts_shorter2_home_female <- Nchelenge5_clusts_shorter2_home[which(Nchelenge5_clusts_shorter2_home$male == 0),] 
Nchelenge5_clusts_shorter2_home_male <- Nchelenge5_clusts_shorter2_home[which(Nchelenge5_clusts_shorter2_home$male == 1),] 

no_home_female <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_home_female$partid))
## one female with no home
no_home_male <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_home_male$partid))
## no males with no home

### number of 'trips' to home loc ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2_home3_female <- Nchelenge5_clust_trips_shorter2_home3[which(Nchelenge5_clust_trips_shorter2_home3$male == 0),]
Nchelenge5_clust_trips_shorter2_home3_male <- Nchelenge5_clust_trips_shorter2_home3[which(Nchelenge5_clust_trips_shorter2_home3$male == 1),]

favstats(Nchelenge5_clust_trips_shorter2_home3_female$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#   1  3      7 13  33 10.08889 9.206772 45       0
favstats(Nchelenge5_clust_trips_shorter2_home3_male$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#      1  4      9 14.75  48 11.13333 10.25783 30       0

### total time at home loc ####
Nchelenge5_clust_trips_shorter2_home4_female <- Nchelenge5_clust_trips_shorter2_home4[which(Nchelenge5_clust_trips_shorter2_home4$male == 0),]
Nchelenge5_clust_trips_shorter2_home4_male <- Nchelenge5_clust_trips_shorter2_home4[which(Nchelenge5_clust_trips_shorter2_home4$male == 1),]

favstats(Nchelenge5_clust_trips_shorter2_home4_female$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
# 0.099914 17.22422 24.67409 29.96997 40.73759 23.06398 9.917224 45       0
favstats(Nchelenge5_clust_trips_shorter2_home4_male$clust_trips_time_new)
#         min       Q1   median       Q3      max     mean       sd  n missing
# 0.028894 12.3519 21.30148 26.13971 36.04353 19.72536 9.983579 30       0
#
### percent time at home ####
favstats(Nchelenge5_clusts_shorter2_home_female$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
#  11.62073 82.43695 92.01757 97.95048 100 84.76494 21.34783 44       0
favstats(Nchelenge5_clusts_shorter2_home_male$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#  1.618541 80.09204 89.52041 98.0249 100 83.26752 20.69486 30       0

perc_time_home2_gender <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_home_female$perc_time_home, Nchelenge5_clusts_shorter2_home_male$perc_time_home),
                                              c(rep("female",length(Nchelenge5_clusts_shorter2_home_female$perc_time_home)),
                                                rep("male",length(Nchelenge5_clusts_shorter2_home_male$perc_time_home)))))
colnames(perc_time_home2_gender) <- c("values","type")
perc_time_home2_gender$values <- as.numeric(as.character(perc_time_home2_gender$values))
perc_time_home2_gender$type <- factor(perc_time_home2_gender$type)

kruskal.test(perc_time_home2_gender$values, perc_time_home2_gender$type) # p = 0.53

ggplot(data=perc_time_home2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,102))
##
### Average Direction of Travel #####
Nchelenge5_clusts_shorter2_no_home$male <- as.factor(Nchelenge5_clusts_shorter2_no_home$male)
temp <- Nchelenge5_clusts_shorter2_no_home[!duplicated(Nchelenge5_clusts_shorter2_no_home$partid),]
temp2 <- temp[which(!(is.na(temp$male))),]

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=male), position=position_dodge2(preserve = "single"), width=4) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Male", values=c("deeppink","blue"))
#

favstats(temp2$avg_dir_from_home_bearing~temp2$male)
#   temp2$male       min       Q1   median       Q3      max     mean       sd  n missing
# 1          0  7.761391 122.3533 208.6944 294.8238 358.9609 206.2246 95.68891 41       0
# 2          1 15.655534 135.8979 223.5365 290.7518 339.2409 208.1619 94.39310 28       0
favstats(temp2$avg_dir_from_home_mag~temp2$male)
#   temp2$male       min        Q1   median       Q3      max     mean       sd  n missing
# 1          0 0.1811803 0.6055279 1.190167 3.681733 21.56929 3.638021 5.447389 41       0
# 2          1 0.1034942 0.5064976 1.662227 3.177905 19.67457 3.448438 4.983412 28       0

temp3 <- temp2[which(temp2$avg_dir_from_home_mag < 65),] ## 164 of 176 parts

favstats(temp3$avg_dir_from_home_bearing~temp3$male)
#   temp3$male       min       Q1   median       Q3      max     mean       sd  n missing
# 1          0  7.761391 122.3533 208.6944 294.8238 358.9609 206.2246 95.68891 41       0
# 2          1 15.655534 135.8979 223.5365 290.7518 339.2409 208.1619 94.39310 28       0
favstats(temp3$avg_dir_from_home_mag~temp3$male)
#   temp3$male       min        Q1   median       Q3      max     mean       sd  n missing
# 1          0 0.1811803 0.6055279 1.190167 3.681733 21.56929 3.638021 5.447389 41       0
# 2          1 0.1034942 0.5064976 1.662227 3.177905 19.67457 3.448438 4.983412 28       0
kruskal.test(temp3$avg_dir_from_home_bearing, temp3$male) # p = 0.53

ggplot(data=temp3, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=male), position=position_dodge2(preserve = "single"), width=4) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Male", values=c("red","blue"))


temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)] <- -(360-temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)])
ggplot(data=temp2, aes(x=avg_dir_from_home_values, y=..count.., fill=male)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~male, labeller = labeller(male = c("0"="Female", "1"="Male")))+
  scale_fill_manual(values=c("red","blue")) +
  scale_y_continuous(breaks=c(seq(0,14,2)))

#

### time spent at clusters vs elsewhere (travel, less important places, etc) ####
Nchelenge5_clusts_shorter2_female <- Nchelenge5_clusts_shorter2[which(Nchelenge5_clusts_shorter2$male == 0),]
Nchelenge5_clusts_shorter2_male <- Nchelenge5_clusts_shorter2[which(Nchelenge5_clusts_shorter2$male == 1),]
Nchelenge5_clusts_shorter2_female$perc_time_in_all_clusts <- Nchelenge5_clusts_shorter2_female$part_trips_time_new/Nchelenge5_clusts_shorter2_female$total_time
Nchelenge5_clusts_shorter2_male$perc_time_in_all_clusts <- Nchelenge5_clusts_shorter2_male$part_trips_time_new/Nchelenge5_clusts_shorter2_male$total_time

summary(Nchelenge5_clusts_shorter2_female$perc_time_in_all_clusts[!duplicated(Nchelenge5_clusts_shorter2_female$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7556  0.8992  0.9560  0.9374  0.9872  0.9984 
summary(Nchelenge5_clusts_shorter2_male$perc_time_in_all_clusts[!duplicated(Nchelenge5_clusts_shorter2_male$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6629  0.8621  0.9165  0.8969  0.9323  0.9955 

# time outside clusters
summary(Nchelenge5_clusts_shorter2_female$total_time[!duplicated(Nchelenge5_clusts_shorter2_female$partid)] - Nchelenge5_clusts_shorter2_female$part_trips_time_new[!duplicated(Nchelenge5_clusts_shorter2_female$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1582  0.5242  1.3370  2.1695  3.8729  8.3344 
summary(Nchelenge5_clusts_shorter2_male$total_time[!duplicated(Nchelenge5_clusts_shorter2_male$partid)] - Nchelenge5_clusts_shorter2_male$part_trips_time_new[!duplicated(Nchelenge5_clusts_shorter2_male$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1395  1.6176  2.8666  3.3659  4.7557  9.8070 
#
### Time at locations ####
favstats(Nchelenge5_clusts_shorter2_no_home_female$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2502778 0.9848611 2.277083 5.764653 784.0228 16.63359 69.54821 368       0
favstats(Nchelenge5_clusts_shorter2_no_home_male$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2425 1.139583 2.211111 4.315139 448.2244 10.93217 37.32198 275       0

### PER PART ###
clust_time_part_median_female <- (favstats(Nchelenge5_clusts_shorter2_no_home_female$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_female$partid)$median)*24
clust_time_part_median_male <- (favstats(Nchelenge5_clusts_shorter2_no_home_male$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_male$partid)$median)*24
favstats(clust_time_part_median_female)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4605556 1.655243 2.381389 3.220312 12.29667 2.827497 2.054988 42       0
favstats(clust_time_part_median_male)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.8686111 1.549722 2.0575 3.195833 17.23194 3.175048 3.475137 29       0

clust_time_part_gender <- as.data.frame(cbind(c(clust_time_part_median_female, clust_time_part_median_male),
                                              c(rep("female",length(clust_time_part_median_female)),
                                                rep("male",length(clust_time_part_median_male)))))
colnames(clust_time_part_gender) <- c("values","type")
clust_time_part_gender$values <- as.numeric(as.character(clust_time_part_gender$values))
clust_time_part_gender$type <- factor(clust_time_part_gender$type)

kruskal.test(clust_time_part_gender$values, clust_time_part_gender$type) # p = 0.76

ggplot(data=clust_time_part_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-1,20)) +
  scale_x_continuous(breaks=c(seq(0,20,2)))



clust_time_part_gender <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_no_home_female$clust_trips_time_new*24, Nchelenge5_clusts_shorter2_no_home_male$clust_trips_time_new*24),
                                        c(rep("female",length(Nchelenge5_clusts_shorter2_no_home_female$clust_trips_time_new)),
                                          rep("male",length(Nchelenge5_clusts_shorter2_no_home_male$clust_trips_time_new)))))
colnames(clust_time_part_gender) <- c("values","type")
clust_time_part_gender$values <- as.numeric(as.character(clust_time_part_gender$values))
clust_time_part_gender$type <- factor(clust_time_part_gender$type)

kruskal.test(clust_time_part_gender$values, clust_time_part_gender$type) # p = 0.99

ggplot(data=clust_time_part_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-1,20)) +
  scale_x_continuous(breaks=c(seq(0,20,2)))

#
#
#
#
####### Biting time metrics #####
### number locations outside home ####
Nchelenge5_clusts_shorter2_no_home_biting_places2_female <- Nchelenge5_clusts_shorter2_no_home_biting_places2[which(Nchelenge5_clusts_shorter2_no_home_biting_places2$male == 0),]
Nchelenge5_clusts_shorter2_no_home_biting_places2_male <- Nchelenge5_clusts_shorter2_no_home_biting_places2[which(Nchelenge5_clusts_shorter2_no_home_biting_places2$male == 1),]

loc_counts_female <- favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_female$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_female$partid)$n
loc_counts_male <- favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_male$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_male$partid)$n
only_home_female <- nlevels(as.factor(Nchelenge5_clusts_shorter2_female$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_female$partid))
only_home_male <- nlevels(as.factor(Nchelenge5_clusts_shorter2_male$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_male$partid))
loc_counts_female2 <- c(loc_counts_female,rep(0, only_home_female))
loc_counts_male2 <- c(loc_counts_male,rep(0, only_home_male))

favstats(loc_counts_female2)
# min Q1 median Q3 max     mean       sd  n missing
#   0  0      1  2   8 1.511 1.996 45       0
favstats(loc_counts_male2)
# min Q1 median Q3 max     mean       sd  n missing
#   0  0    1.5 2.75  12 2.033 2.697 30       0
summary(as.factor(loc_counts_female2))
#  0  1  2  3  4  6  8 
# 18 13  5  1  4  3  1 
summary(as.factor(loc_counts_male2))
# 0  1  2  3  4  5  9 12 
# 10  5  7  4  1  1  1  1 
## all those with 0 have biting time at home location (not just nowhere)

locs_count2 <- as.data.frame(cbind(c(loc_counts_female2, loc_counts_male2),
                                   c(rep("female",length(loc_counts_female2)),
                                     rep("male",length(loc_counts_male2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-0.35,12.35)) +
  scale_x_continuous(breaks=c(seq(0,12,1)))
#
##
### time per location ####
favstats((Nchelenge5_clusts_shorter2_no_home_biting_places2_female$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5112 1.693  7.973 26.2 285.4 25.35 47.59 68       0
favstats((Nchelenge5_clusts_shorter2_no_home_biting_places2_male$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5511 1.578  5.434 14.54 193.3 15.41 29.5 61       0

clusts_median_Nchelenge_female <- (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_female$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_female$partid)$median)*24
clusts_median_Nchelenge_male <- (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_male$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_male$partid)$median)*24
favstats(clusts_median_Nchelenge_female)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.8618 3.712  11.04 20.76 79.21 16.03 17.49 27       0
favstats(clusts_median_Nchelenge_male)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.8406 2.193  8.385 27.76 82.97 17.91 22.85 20       0

clust_time_part_gender <- as.data.frame(cbind(c(clusts_median_Nchelenge_female, clusts_median_Nchelenge_male),
                                              c(rep("female",length(clusts_median_Nchelenge_female)),
                                                rep("male",length(clusts_median_Nchelenge_male)))))
colnames(clust_time_part_gender) <- c("values","type")
clust_time_part_gender$values <- as.numeric(as.character(clust_time_part_gender$values))
clust_time_part_gender$type <- factor(clust_time_part_gender$type)

ggplot(data=clust_time_part_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-2,100)) +
  scale_x_continuous(breaks=c(seq(0,100,20)))
#
### number trips ####
Nchelenge5_clust_trips_shorter2_no_home_biting_places3_female <- Nchelenge5_clust_trips_shorter2_no_home_biting_places3[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$male == 0),]
Nchelenge5_clust_trips_shorter2_no_home_biting_places3_male <- Nchelenge5_clust_trips_shorter2_no_home_biting_places3[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$male == 1),]

favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2  10 1.985 1.973 67       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2  13 2.237 2.452 59       0

summary(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2))
# 1  2  3  4  5  7  9 10 
# 43 10  6  3  1  1  2  1 
summary(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2))
#  1  2  3  5  6  7 12 13 
# 34 12  6  1  2  2  1  1 

### PER PART ###
trip_count_median_female <- favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_female$partid)$median
favstats(trip_count_median_female)
# min Q1 median Q3 max  mean  sd  n missing
#    1  1      1 1.75   7 1.63 1.305 27       0
trip_count_median_male <- favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_male$partid)$median
favstats(trip_count_median_male)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1    1.5  3 4.5 1.95 1.099 20       0

trips_part2_gender <- as.data.frame(cbind(c(trip_count_median_female, trip_count_median_male),
                                          c(rep("female",length(trip_count_median_female)),
                                            rep("male",length(trip_count_median_male)))))
colnames(trips_part2_gender) <- c("values","type")
trips_part2_gender$values <- as.numeric(as.character(trips_part2_gender$values))
trips_part2_gender$type <- factor(trips_part2_gender$type)

ggplot(data=trips_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.3) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,8.25)) + 
  scale_x_continuous(breaks=c(seq(0,10,2)))
#
##### percent biting time spent at home vs elsewhere #####
Nchelenge5_clusts_shorter2_biting_places3b_female <- Nchelenge5_clusts_shorter2_biting_places3b[which(Nchelenge5_clusts_shorter2_biting_places3b$male == 0),]
Nchelenge5_clusts_shorter2_biting_places3b_male <- Nchelenge5_clusts_shorter2_biting_places3b[which(Nchelenge5_clusts_shorter2_biting_places3b$male == 1),]

favstats(Nchelenge5_clusts_shorter2_biting_places3b_female$percent_home_biting_time)
favstats(Nchelenge5_clusts_shorter2_biting_places3b_female$percent_home_biting_time[which(Nchelenge5_clusts_shorter2_biting_places3b_female$percent_home_biting_time != 0)])
favstats(Nchelenge5_clusts_shorter2_biting_places3b_male$percent_home_biting_time)
favstats(Nchelenge5_clusts_shorter2_biting_places3b_male$percent_home_biting_time[which(Nchelenge5_clusts_shorter2_biting_places3b_male$percent_home_biting_time != 0)])


perc_time_home2_gender <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_biting_places3b_female$percent_home_biting_time, Nchelenge5_clusts_shorter2_biting_places3b_male$percent_home_biting_time),
                                              c(rep("female",length(Nchelenge5_clusts_shorter2_biting_places3b_female$percent_home_biting_time)),
                                                rep("male",length(Nchelenge5_clusts_shorter2_biting_places3b_male$percent_home_biting_time)))))
colnames(perc_time_home2_gender) <- c("values","type")
perc_time_home2_gender$values <- as.numeric(as.character(perc_time_home2_gender$values))
perc_time_home2_gender$type <- factor(perc_time_home2_gender$type)

ggplot(data=perc_time_home2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,100))
##
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
###### ###### AGE SPLIT ###### ######
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
###############
###############
##### CHOMA #####
###### General Age Analysis ####
Choma_part_info <- Choma5_clusts_shorter2[!duplicated(Choma5_clusts_shorter2$partid),c(1:3)]
summary(Choma_part_info$age_all)
# 8 - 78; mean=38, median=39, 1Q=20, 3Q=51.8

Choma_part_info$partid2 <- gsub("b","", Choma_part_info$partid)
Choma_part_info2 <- Choma_part_info[!duplicated(Choma_part_info$partid2),]
favstats(as.numeric(Choma_part_info2$age_all))


Choma_part_info$age_18 <- factor(0, levels=c(0,1))
Choma_part_info$age_18[which(Choma_part_info$age_all >= 18)] <- 1
summary(as.factor(Choma_part_info$age_18))
#  0   1
# 12  50  
Choma_part_info$age_16_35_55_up <- factor(0, levels=c(0,1,2,3)) # 0 = <16; 1 = 17-34; 2 = 35-54; 3= 55+
Choma_part_info$age_16_35_55_up[which(Choma_part_info$age_all >= 17 & Choma_part_info$age_all < 35)] <- 1
Choma_part_info$age_16_35_55_up[which(Choma_part_info$age_all >= 35 & Choma_part_info$age_all < 55)] <- 2
Choma_part_info$age_16_35_55_up[which(Choma_part_info$age_all >= 55)] <- 3
summary(as.factor(Choma_part_info$age_16_35_55_up))
#  0   1   2   3
# 10  16  22  14

ggplot(data=Choma_part_info, aes(x=age_all)) +
  geom_histogram(binwidth = 1, aes(fill=age_16_35_55_up)) +
  coord_cartesian(xlim=c(0,80))
#
###### Split by Age ######
Choma_parts_1_16 <- Choma_part_info$partid[which(Choma_part_info$age_16_35_55_up == 0)]
Choma_parts_17_34 <- Choma_part_info$partid[which(Choma_part_info$age_16_35_55_up == 1)]
Choma_parts_35_54 <- Choma_part_info$partid[which(Choma_part_info$age_16_35_55_up == 2)]
Choma_parts_55_up <- Choma_part_info$partid[which(Choma_part_info$age_16_35_55_up == 3)]
##
### number of locations (without home) ####
Choma5_clusts_shorter2_no_home_1_16 <- Choma5_clusts_shorter2_no_home[which(Choma5_clusts_shorter2_no_home$partid %in% Choma_parts_1_16),]
Choma5_clusts_shorter2_no_home_17_34 <- Choma5_clusts_shorter2_no_home[which(Choma5_clusts_shorter2_no_home$partid %in% Choma_parts_17_34),]
Choma5_clusts_shorter2_no_home_35_54 <- Choma5_clusts_shorter2_no_home[which(Choma5_clusts_shorter2_no_home$partid %in% Choma_parts_35_54),]
Choma5_clusts_shorter2_no_home_55_up <- Choma5_clusts_shorter2_no_home[which(Choma5_clusts_shorter2_no_home$partid %in% Choma_parts_55_up),]

loc_counts_1_16 <- favstats(Choma5_clusts_shorter2_no_home_1_16$clust~Choma5_clusts_shorter2_no_home_1_16$partid)$n
only_home_1_16 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_1_16)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_1_16$partid))
loc_counts_1_162 <- c(loc_counts_1_16,rep(0, only_home_1_16))
favstats(loc_counts_1_162)
# min   Q1 median     Q3  max     mean        sd  n missing
#   6  8     13 17  18 12.54545 4.655398 11       0
loc_counts_17_34 <- favstats(Choma5_clusts_shorter2_no_home_17_34$clust~Choma5_clusts_shorter2_no_home_17_34$partid)$n
only_home_17_34 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_17_34)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_17_34$partid))
loc_counts_17_342 <- c(loc_counts_17_34,rep(0, only_home_17_34))
favstats(loc_counts_17_342)
# min   Q1 median     Q3  max     mean        sd  n missing
#   0  5      7 17  20 10.26667 6.902036 15       0
loc_counts_35_54 <- favstats(Choma5_clusts_shorter2_no_home_35_54$clust~Choma5_clusts_shorter2_no_home_35_54$partid)$n
only_home_35_54 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_35_54)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_35_54$partid))
loc_counts_35_542 <- c(loc_counts_35_54,rep(0, only_home_35_54))
favstats(loc_counts_35_542)
# min   Q1 median     Q3  max     mean        sd  n missing
#  3 9.25     12 16.75  22 12.63636 5.778 22       0
loc_counts_55_up <- favstats(Choma5_clusts_shorter2_no_home_55_up$clust~Choma5_clusts_shorter2_no_home_55_up$partid)$n
only_home_55_up <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_55_up)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_55_up$partid))
loc_counts_55_up2 <- c(loc_counts_55_up,rep(0, only_home_55_up))
favstats(loc_counts_55_up2)
# min   Q1 median     Q3  max     mean        sd  n missing
#   5 8.25     12 16  21   12 5.276946 14       0

agecols <- brewer.pal(5, "YlGn")

locs_count2 <- as.data.frame(cbind(c(loc_counts_1_162, loc_counts_17_342, loc_counts_35_542, loc_counts_55_up2),
                                   c(rep("1-16",length(loc_counts_1_162)),
                                     rep("17-34",length(loc_counts_17_342)),
                                     rep("35-54",length(loc_counts_35_542)),
                                     rep("55+",length(loc_counts_55_up2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

shapiro.test(locs_count2$values) ## not normal
kruskal.test(locs_count2$values, locs_count2$type) # p = 0.69
#

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-1,25)) +
  scale_x_continuous(breaks=c(seq(0,25,5)))
#

### distance of locations (without home) ####
### TOTAL ###
favstats(Choma5_clusts_shorter2_no_home_1_16$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
# 0.1001338 0.7454768 1.774908 5.803411 74.90588 5.958273 10.85805 138       0
favstats(Choma5_clusts_shorter2_no_home_17_34$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 0.06883561 1.212972 11.39114 16.88285 79.10662 12.55057 15.10268 154       0
favstats(Choma5_clusts_shorter2_no_home_35_54$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#  0.0343186 1.488736 10.99185 16.78871 112.2062 12.14961 15.36231 278       0
favstats(Choma5_clusts_shorter2_no_home_55_up$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 0.03054886 1.220107 11.65792 16.99273 46.75496 11.04315 9.602567 168       0

### PER PART ###
medians_km_1_16 <- (favstats(Choma5_clusts_shorter2_no_home_1_16$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_1_16$partid)$median)
medians_km_17_34 <- (favstats(Choma5_clusts_shorter2_no_home_17_34$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_17_34$partid)$median)
medians_km_35_54 <- (favstats(Choma5_clusts_shorter2_no_home_35_54$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_35_54$partid)$median)
medians_km_55_up <- (favstats(Choma5_clusts_shorter2_no_home_55_up$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_55_up$partid)$median)

favstats(medians_km_1_16) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
# 0.6884276 0.8466509 1.128486 4.185754 22.62318 4.452841 6.770917 11       0
favstats(medians_km_17_34) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.4184579 0.914198 1.796451 13.64822 18.21496 6.407055 6.825852 14       0
favstats(medians_km_35_54) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.3318577 1.259387 7.721501 14.4196 32.16992 9.532807 8.779758 22       0
favstats(medians_km_55_up) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.2448758 1.399819 13.21831 16.23204 22.49014 10.44794 7.774949 14       0

dists_km2_age <- as.data.frame(cbind(c(medians_km_1_16, medians_km_17_34, medians_km_35_54, medians_km_55_up),
                                     c(rep("1-16",length(medians_km_1_16)),
                                       rep("17-34",length(medians_km_17_34)),
                                       rep("35-54",length(medians_km_35_54)),
                                       rep("55+",length(medians_km_55_up)))))
colnames(dists_km2_age) <- c("values","type")
dists_km2_age$values <- as.numeric(as.character(dists_km2_age$values))
dists_km2_age$type <- factor(dists_km2_age$type)

shapiro.test(dists_km2_age$values) ## not normal
kruskal.test(dists_km2_age$values, dists_km2_age$type) # p = 0.23
pairwise.wilcox.test(dists_km2_age$values, dists_km2_age$type, p.adjust.method = "holm")

ggplot(data=dists_km2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-2,33))
###
### number of trips (without home) ####
### TOTAL ###
Choma5_clust_trips_shorter2_no_home3_1_16 <- Choma5_clust_trips_shorter2_no_home3[which(Choma5_clust_trips_shorter2_no_home3$partid %in% Choma_parts_1_16),]
Choma5_clust_trips_shorter2_no_home3_17_34 <- Choma5_clust_trips_shorter2_no_home3[which(Choma5_clust_trips_shorter2_no_home3$partid %in% Choma_parts_17_34),]
Choma5_clust_trips_shorter2_no_home3_35_54 <- Choma5_clust_trips_shorter2_no_home3[which(Choma5_clust_trips_shorter2_no_home3$partid %in% Choma_parts_35_54),]
Choma5_clust_trips_shorter2_no_home3_55_up <- Choma5_clust_trips_shorter2_no_home3[which(Choma5_clust_trips_shorter2_no_home3$partid %in% Choma_parts_55_up),]

favstats(Choma5_clust_trips_shorter2_no_home3_1_16$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      2  3  35 3.224638 5.177179 138       0
favstats(Choma5_clust_trips_shorter2_no_home3_17_34$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
# 1  1      1  2  30 2.616883 4.238322 154       0
favstats(Choma5_clust_trips_shorter2_no_home3_35_54$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      2  3  52 3.068345 5.045533 278       0
favstats(Choma5_clust_trips_shorter2_no_home3_55_up$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      2  3  58 3.333333 6.487821 168       0

### PER PART ###
trips_part_median_1_16 <- favstats(Choma5_clust_trips_shorter2_no_home3_1_16$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_1_16$partid)$median
trips_part_median_17_34 <- favstats(Choma5_clust_trips_shorter2_no_home3_17_34$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_17_34$partid)$median
trips_part_median_35_54 <- favstats(Choma5_clust_trips_shorter2_no_home3_35_54$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_35_54$partid)$median
trips_part_median_55_up <- favstats(Choma5_clust_trips_shorter2_no_home3_55_up$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_55_up$partid)$median

favstats(trips_part_median_1_16)
#   min  Q1   median   Q3      max       mean       sd  n missing
#    1  1      2  2   4 1.727273 0.904534 11       0
favstats(trips_part_median_17_34)
#   min Q1 median  Q3   max      mean       sd  n missing
#   1  1      1 1.375   3 1.321429 0.6078696 14       0
favstats(trips_part_median_35_54)
#   min  Q1   median   Q3      max       mean       sd  n missing
#   1  1      1  2   3 1.477273 0.587109 22       0
favstats(trips_part_median_55_up)
#   min Q1 median  Q3   max      mean       sd  n missing
#   1  1      2  2   3 1.714286 0.6112498 14       0

trips_part2_age <- as.data.frame(cbind(c(trips_part_median_1_16, trips_part_median_17_34, trips_part_median_35_54, trips_part_median_55_up),
                                       c(rep("1-16",length(trips_part_median_1_16)),
                                         rep("17-34",length(trips_part_median_17_34)),
                                         rep("35-54",length(trips_part_median_35_54)),
                                         rep("55+",length(trips_part_median_55_up)))))
colnames(trips_part2_age) <- c("values","type")
trips_part2_age$values <- as.numeric(as.character(trips_part2_age$values))
trips_part2_age$type <- factor(trips_part2_age$type)
kruskal.test(trips_part2_age$values, trips_part2_age$type) # p = 0.16
pairwise.wilcox.test(trips_part2_age$values, trips_part2_age$type, p.adjust.method = "holm")

ggplot(data=trips_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2.5) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,6)) + 
  scale_x_continuous(breaks=c(seq(0,6,2)))
#
### time per trip ####
### TOTAL ###
Choma5_clust_trips_shorter2_no_home_1_16 <- Choma5_clust_trips_shorter2_no_home[which(Choma5_clust_trips_shorter2_no_home$partid %in% Choma_parts_1_16),]
Choma5_clust_trips_shorter2_no_home_17_34 <- Choma5_clust_trips_shorter2_no_home[which(Choma5_clust_trips_shorter2_no_home$partid %in% Choma_parts_17_34),]
Choma5_clust_trips_shorter2_no_home_35_54 <- Choma5_clust_trips_shorter2_no_home[which(Choma5_clust_trips_shorter2_no_home$partid %in% Choma_parts_35_54),]
Choma5_clust_trips_shorter2_no_home_55_up <- Choma5_clust_trips_shorter2_no_home[which(Choma5_clust_trips_shorter2_no_home$partid %in% Choma_parts_55_up),]

favstats(Choma5_clust_trips_shorter2_no_home_1_16$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2408333 0.7392361   1.44 3.359931 472.56 5.907219 25.81055 434       0
favstats(Choma5_clust_trips_shorter2_no_home_17_34$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2402778 0.6031944 1.197778 3.684167 206.0228 8.361131 25.38287 399       0
favstats(Choma5_clust_trips_shorter2_no_home_35_54$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2402778 0.6679861 1.546389 4.313889 355.0831 8.229321 25.30912 848       0
favstats(Choma5_clust_trips_shorter2_no_home_55_up$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2491667 0.6856944   1.44 4.804306 456.055 13.17145 46.8489 555       0

### PER PART ###
trip_time_part_median_1_16 <- favstats(Choma5_clust_trips_shorter2_no_home_1_16$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_1_16$partid)$median
trip_time_part_median_17_34 <- favstats(Choma5_clust_trips_shorter2_no_home_17_34$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_17_34$partid)$median
trip_time_part_median_35_54 <- favstats(Choma5_clust_trips_shorter2_no_home_35_54$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_35_54$partid)$median
trip_time_part_median_55_up <- favstats(Choma5_clust_trips_shorter2_no_home_55_up$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_55_up$partid)$median

favstats(trip_time_part_median_1_16)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.7469444 1.208611 1.319861 1.547361 3.36 1.495177 0.6881796 11       0
favstats(trip_time_part_median_17_34)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.5879167 0.8911111 1.186667 1.489201 2.1575 1.215446 0.4487393 14       0
favstats(trip_time_part_median_35_54)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.7531944 1.275069 1.4775 1.734097 4.093333 1.604861 0.7363699 22       0
favstats(trip_time_part_median_55_up)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.6075 1.354062 1.584653 2.462917 5.650556 2.219216 1.495475 14       0

trip_time_part_median_55_up_short <- trip_time_part_median_55_up ## no outliers to remove

trip_times_part2_age <- as.data.frame(cbind(c(trip_time_part_median_1_16, trip_time_part_median_17_34, trip_time_part_median_35_54, trip_time_part_median_55_up_short),
                                            c(rep("1-16",length(trip_time_part_median_1_16)),
                                              rep("17-34",length(trip_time_part_median_17_34)),
                                              rep("35-54",length(trip_time_part_median_35_54)),
                                              rep("55+",length(trip_time_part_median_55_up_short)))))
colnames(trip_times_part2_age) <- c("values","type")
trip_times_part2_age$values <- as.numeric(as.character(trip_times_part2_age$values))
trip_times_part2_age$type <- factor(trip_times_part2_age$type)
kruskal.test(trip_times_part2_age$values, trip_times_part2_age$type) # p = 0.06
pairwise.wilcox.test(trip_times_part2_age$values, trip_times_part2_age$type, p.adjust.method = "holm")

ggplot(data=trip_times_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,8))

#
### Home Loc ####
Choma5_clusts_shorter2_home_1_16 <- Choma5_clusts_shorter2_home[which(Choma5_clusts_shorter2_home$partid %in% Choma_parts_1_16),] 
no_home_1_16 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_1_16)])) - nlevels(as.factor(Choma5_clusts_shorter2_home_1_16$partid))# 1
## 2 with no home
Choma5_clusts_shorter2_home_17_34 <- Choma5_clusts_shorter2_home[which(Choma5_clusts_shorter2_home$partid %in% Choma_parts_17_34),] 
no_home_17_34 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_17_34)])) - nlevels(as.factor(Choma5_clusts_shorter2_home_17_34$partid))# 1
## none with no home
Choma5_clusts_shorter2_home_35_54 <- Choma5_clusts_shorter2_home[which(Choma5_clusts_shorter2_home$partid %in% Choma_parts_35_54),] 
no_home_35_54 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_35_54)])) - nlevels(as.factor(Choma5_clusts_shorter2_home_35_54$partid))# 1
## none with no home
Choma5_clusts_shorter2_home_55_up <- Choma5_clusts_shorter2_home[which(Choma5_clusts_shorter2_home$partid %in% Choma_parts_55_up),] 
no_home_55_up <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_55_up)])) - nlevels(as.factor(Choma5_clusts_shorter2_home_55_up$partid))# 1
## 1 with no home

### number of 'trips' to home loc ####
### TOTAL ###
Choma5_clust_trips_shorter2_home3_1_16 <- Choma5_clust_trips_shorter2_home3[which(Choma5_clust_trips_shorter2_home3$partid %in% Choma_parts_1_16),]
Choma5_clust_trips_shorter2_home3_17_34 <- Choma5_clust_trips_shorter2_home3[which(Choma5_clust_trips_shorter2_home3$partid %in% Choma_parts_17_34),]
Choma5_clust_trips_shorter2_home3_35_54 <- Choma5_clust_trips_shorter2_home3[which(Choma5_clust_trips_shorter2_home3$partid %in% Choma_parts_35_54),]
Choma5_clust_trips_shorter2_home3_55_up <- Choma5_clust_trips_shorter2_home3[which(Choma5_clust_trips_shorter2_home3$partid %in% Choma_parts_55_up),]

favstats(Choma5_clust_trips_shorter2_home3_1_16$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#   2  5     15 26  37   17 12.98075 9       0
favstats(Choma5_clust_trips_shorter2_home3_17_34$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#    1  4      6 10.5  20  7.8 5.294202 15       0
favstats(Choma5_clust_trips_shorter2_home3_35_54$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#   1  4      6 12  38 9.954545 10.04503 22       0
favstats(Choma5_clust_trips_shorter2_home3_55_up$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#  2  3      4 20  37 12.92308 13.81703 13       0
#
### total time at home loc ####
Choma5_clust_trips_shorter2_home4_1_16 <- Choma5_clust_trips_shorter2_home4[which(Choma5_clust_trips_shorter2_home4$partid %in% Choma_parts_1_16),]
favstats(Choma5_clust_trips_shorter2_home4_1_16$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
#  0.1946354 21.04891 24.2914 25.68119 36.24018 20.8416 11.56993 9       0
Choma5_clust_trips_shorter2_home4_17_34 <- Choma5_clust_trips_shorter2_home4[which(Choma5_clust_trips_shorter2_home4$partid %in% Choma_parts_17_34),]
favstats(Choma5_clust_trips_shorter2_home4_17_34$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
#  1.165226 16.31526 23.96804 29.64717 40.365 21.53748 11.71588 15       0
Choma5_clust_trips_shorter2_home4_35_54 <- Choma5_clust_trips_shorter2_home4[which(Choma5_clust_trips_shorter2_home4$partid %in% Choma_parts_35_54),]
favstats(Choma5_clust_trips_shorter2_home4_35_54$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
# 0.02523149 5.496441 18.47617 25.46826 116.9807 21.85043 23.91561 22       0
Choma5_clust_trips_shorter2_home4_55_up <- Choma5_clust_trips_shorter2_home4[which(Choma5_clust_trips_shorter2_home4$partid %in% Choma_parts_55_up),]
favstats(Choma5_clust_trips_shorter2_home4_55_up$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
#  0.0630382 2.006731 7.155851 22.6205 32.67663 12.75947 12.1248 13       0
#
### percent time at home ####
## without 3 people who spent 0% of time at home
favstats(Choma5_clusts_shorter2_home_1_16$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
#  0.7856377 80.73163 89.01833 96.04384 98.82904 73.26005 37.57168 9       0
favstats(Choma5_clusts_shorter2_home_17_34$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#  3.773483 58.35464 91.52463 97.59479 100 72.98601 37.36857 15       0
favstats(Choma5_clusts_shorter2_home_35_54$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
#  0.08520987 21.72893 72.51377 96.1336 99.05443 59.95318 36.93578 22       0
favstats(Choma5_clusts_shorter2_home_55_up$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#    0.243299 2.778517 19.4021 75.7409 94.76292 35.76407 37.80257 13       0

perc_home2_age <- as.data.frame(cbind(c(Choma5_clusts_shorter2_home_1_16$perc_time_home, Choma5_clusts_shorter2_home_17_34$perc_time_home, Choma5_clusts_shorter2_home_35_54$perc_time_home, Choma5_clusts_shorter2_home_55_up$perc_time_home),
                                      c(rep("1-16",length(Choma5_clusts_shorter2_home_1_16$perc_time_home)),
                                        rep("17-34",length(Choma5_clusts_shorter2_home_17_34$perc_time_home)),
                                        rep("35-54",length(Choma5_clusts_shorter2_home_35_54$perc_time_home)),
                                        rep("55+",length(Choma5_clusts_shorter2_home_55_up$perc_time_home)))))
colnames(perc_home2_age) <- c("values","type")
perc_home2_age$values <- as.numeric(as.character(perc_home2_age$values))
perc_home2_age$type <- factor(perc_home2_age$type)
kruskal.test(perc_home2_age$values, perc_home2_age$type) # p = 0.03
pairwise.wilcox.test(perc_home2_age$values, perc_home2_age$type, p.adjust.method = "holm")

ggplot(data=perc_home2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-1,100))
##
### Average Direction of Travel #####
temp <- Choma5_clusts_shorter2_no_home[!duplicated(Choma5_clusts_shorter2_no_home$partid),]

age_group <- Choma5_clusts_shorter2_no_home$age_all[!duplicated(Choma5_clusts_shorter2_no_home$partid)]

temp$age_group <- factor(0, levels=c(0:3)) # 0 = <16; 1 = 17-34; 2 = 35-54; 3= 55+
temp$age_group[which(age_group >= 17 & age_group < 35)] <- 1
temp$age_group[which(age_group >= 35 & age_group < 55)] <- 2
temp$age_group[which(age_group >= 55)] <- 3
summary(as.factor(temp$age_group))

favstats(temp$avg_dir_from_home_bearing~temp$age_group)
# temp$age_group      min       Q1    median       Q3      max     mean        sd  n missing
# 1              3 21.03243 70.28968  75.70738 146.3676 302.3886 111.7892  79.81307 14       0
# 2              2 21.74463 70.14785 119.25977 232.1584 308.4303 150.2771  94.77611 22       0
# 3              1 16.56593 74.70586 146.87705 280.4911 334.6366 165.5802 117.80319 14       0
# 4              0  8.43748 48.02399 122.31646 272.5941 344.5980 152.4312 120.33806 11       0
favstats(temp$avg_dir_from_home_mag~temp$age_group)
# temp$age_group       min        Q1    median        Q3      max      mean       sd  n
# 1              3 0.5553266 3.9904048 12.872336 13.761859 19.25974  9.902251 6.142657 14       0
# 2              2 0.3610027 1.2553177  9.397531 16.046085 33.03357 10.389801 9.735389 22       0
# 3              1 0.1915475 0.8755871  3.502125 14.757321 22.20797  8.047884 8.506031 14       0
# 4              0 0.2044328 0.6244360  1.663715  6.102087 22.15730  4.932353 7.047860 11       0

kruskal.test(temp$avg_dir_from_home_bearing, temp$age_group) # p = 0.67
pairwise.wilcox.test(temp$avg_dir_from_home_bearing, temp$age_group, p.adjust.method = "holm")


ggplot(data=temp, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=age_group), position=position_dodge2(preserve = "single"), width=4) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Age Group", values=agecols[c(5,4,2,1)], labels=c("55+", "35-54", "17-34", "0-16"))





temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)] <- -(360-temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)])
ggplot(data=temp, aes(x=avg_dir_from_home_values, y=..count.., fill=age_group)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~age_group, labeller = labeller(age_group = c("0"="0-16", "1"="17-34","2"="35-54", "3"="55+")))+
  scale_fill_manual(values=agecols[c(1,2,4,5)]) 

#
### time spent at clusters vs elsewhere (travel, less important places, etc) ####
# percent time in clusters
Choma5_clusts_shorter2_1_16 <- Choma5_clusts_shorter2[which(Choma5_clusts_shorter2$partid %in% Choma_parts_1_16),]
Choma5_clusts_shorter2_17_34 <- Choma5_clusts_shorter2[which(Choma5_clusts_shorter2$partid %in% Choma_parts_17_34),]
Choma5_clusts_shorter2_35_54 <- Choma5_clusts_shorter2[which(Choma5_clusts_shorter2$partid %in% Choma_parts_35_54),]
Choma5_clusts_shorter2_55_up <- Choma5_clusts_shorter2[which(Choma5_clusts_shorter2$partid %in% Choma_parts_55_up),]

Choma5_clusts_shorter2_1_16$perc_time_in_all_clusts <- Choma5_clusts_shorter2_1_16$part_trips_time_new/Choma5_clusts_shorter2_1_16$total_time
Choma5_clusts_shorter2_17_34$perc_time_in_all_clusts <- Choma5_clusts_shorter2_17_34$part_trips_time_new/Choma5_clusts_shorter2_17_34$total_time
Choma5_clusts_shorter2_35_54$perc_time_in_all_clusts <- Choma5_clusts_shorter2_35_54$part_trips_time_new/Choma5_clusts_shorter2_35_54$total_time
Choma5_clusts_shorter2_55_up$perc_time_in_all_clusts <- Choma5_clusts_shorter2_55_up$part_trips_time_new/Choma5_clusts_shorter2_55_up$total_time

summary(Choma5_clusts_shorter2_1_16$perc_time_in_all_clusts[!duplicated(Choma5_clusts_shorter2_1_16$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8451  0.8610  0.9074  0.9118  0.9541  0.9932 
summary(Choma5_clusts_shorter2_17_34$perc_time_in_all_clusts[!duplicated(Choma5_clusts_shorter2_17_34$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8663  0.9148  0.9509  0.9451  0.9895  0.9955 
summary(Choma5_clusts_shorter2_35_54$perc_time_in_all_clusts[!duplicated(Choma5_clusts_shorter2_35_54$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6629  0.8749  0.9248  0.9012  0.9642  0.9984 
summary(Choma5_clusts_shorter2_55_up$perc_time_in_all_clusts[!duplicated(Choma5_clusts_shorter2_55_up$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7781  0.8729  0.9317  0.9195  0.9782  0.9950 


# time outside clusters
summary(Choma5_clusts_shorter2_1_16$total_time[!duplicated(Choma5_clusts_shorter2_1_16$partid)] - Choma5_clusts_shorter2_1_16$part_trips_time_new[!duplicated(Choma5_clusts_shorter2_1_16$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2509  1.2039  2.9537  2.6542  4.3409  4.8027 
summary(Choma5_clusts_shorter2_17_34$total_time[!duplicated(Choma5_clusts_shorter2_17_34$partid)] - Choma5_clusts_shorter2_17_34$part_trips_time_new[!duplicated(Choma5_clusts_shorter2_17_34$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1395  0.4087  1.4091  1.8793  3.0862  4.7086 
summary(Choma5_clusts_shorter2_35_54$total_time[!duplicated(Choma5_clusts_shorter2_35_54$partid)] - Choma5_clusts_shorter2_35_54$part_trips_time_new[!duplicated(Choma5_clusts_shorter2_35_54$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2275  1.1145  2.0146  3.3552  4.5656  9.8070 
summary(Choma5_clusts_shorter2_55_up$total_time[!duplicated(Choma5_clusts_shorter2_55_up$partid)] - Choma5_clusts_shorter2_55_up$part_trips_time_new[!duplicated(Choma5_clusts_shorter2_55_up$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2383  0.7049  2.3483  2.7655  4.5573  7.0787 
#
### Time at locations ####
favstats(Choma5_clusts_shorter2_no_home_1_16$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.3111111 1.2 1.920694 4.799931 488.4 18.60036 77.23378 138       0
favstats(Choma5_clusts_shorter2_no_home_17_34$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2461111 0.8754167 1.559583 3.228889 885.6553 21.67216 108.6303 154       0
favstats(Choma5_clusts_shorter2_no_home_35_54$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2686111 1.119861 2.120556 4.534931 704.1456 25.10826 104.2773 278       0
favstats(Choma5_clusts_shorter2_no_home_55_up$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2644444 1.139375 2.031528 5.64625 834.3897 43.52203 154.2013 168       0

### PER PART ###
clust_time_part_median_1_16 <- (favstats(Choma5_clusts_shorter2_no_home_1_16$clust_trips_time_new~Choma5_clusts_shorter2_no_home_1_16$partid)$median)*24
clust_time_part_median_17_34 <- (favstats(Choma5_clusts_shorter2_no_home_17_34$clust_trips_time_new~Choma5_clusts_shorter2_no_home_17_34$partid)$median)*24
clust_time_part_median_35_54 <- (favstats(Choma5_clusts_shorter2_no_home_35_54$clust_trips_time_new~Choma5_clusts_shorter2_no_home_35_54$partid)$median)*24
clust_time_part_median_55_up <- (favstats(Choma5_clusts_shorter2_no_home_55_up$clust_trips_time_new~Choma5_clusts_shorter2_no_home_55_up$partid)$median)*24

favstats(clust_time_part_median_1_16)
#       min       Q1   median       Q3      max     mean      sd  n missing
#  1.309583 1.372361 1.923611 2.136319 6.960278 2.229571 1.628297 11       0
favstats(clust_time_part_median_17_34)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.7801389 1.264722 1.531111 2.055174 2.410833 1.603393 0.5096754 14       0
favstats(clust_time_part_median_35_54)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.8297222 1.405 1.849444 2.579792 3.737639 2.023479 0.8074491 22       0
favstats(clust_time_part_median_55_up)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.200278 1.731354 2.031528 3.699583 7.662778 2.841647 1.76496 14       0


clust_time_part_age <- as.data.frame(cbind(c(clust_time_part_median_1_16, clust_time_part_median_17_34,clust_time_part_median_35_54,clust_time_part_median_55_up),
                                              c(rep("1-16",length(clust_time_part_median_1_16)),
                                                rep("17-34",length(clust_time_part_median_17_34)),
                                                rep("35-54",length(clust_time_part_median_35_54)),
                                                rep("55+",length(clust_time_part_median_55_up)))))
colnames(clust_time_part_age) <- c("values","type")
clust_time_part_age$values <- as.numeric(as.character(clust_time_part_age$values))
clust_time_part_age$type <- factor(clust_time_part_age$type)
kruskal.test(clust_time_part_age$values, clust_time_part_age$type) # p = 0.06
pairwise.wilcox.test(clust_time_part_age$values, clust_time_part_age$type, p.adjust.method = "holm")

ggplot(data=clust_time_part_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,9)) +
  scale_x_continuous(breaks=c(seq(0,9,1)))
#
#
#
####### Biting time metrics #####
### number locations outside home ####
Choma5_clusts_shorter2_no_home_biting_places2_1_16 <- Choma5_clusts_shorter2_no_home_biting_places2[which(Choma5_clusts_shorter2_no_home_biting_places2$partid %in% Choma_parts_1_16),]
Choma5_clusts_shorter2_no_home_biting_places2_17_34 <- Choma5_clusts_shorter2_no_home_biting_places2[which(Choma5_clusts_shorter2_no_home_biting_places2$partid %in% Choma_parts_17_34),]
Choma5_clusts_shorter2_no_home_biting_places2_35_54 <- Choma5_clusts_shorter2_no_home_biting_places2[which(Choma5_clusts_shorter2_no_home_biting_places2$partid %in% Choma_parts_35_54),]
Choma5_clusts_shorter2_no_home_biting_places2_55_up <- Choma5_clusts_shorter2_no_home_biting_places2[which(Choma5_clusts_shorter2_no_home_biting_places2$partid %in% Choma_parts_55_up),]

loc_counts_1_16 <- favstats(Choma5_clusts_shorter2_no_home_biting_places2_1_16$clust~Choma5_clusts_shorter2_no_home_biting_places2_1_16$partid)$n
loc_counts_17_34 <- favstats(Choma5_clusts_shorter2_no_home_biting_places2_17_34$clust~Choma5_clusts_shorter2_no_home_biting_places2_17_34$partid)$n
loc_counts_35_54 <- favstats(Choma5_clusts_shorter2_no_home_biting_places2_35_54$clust~Choma5_clusts_shorter2_no_home_biting_places2_35_54$partid)$n
loc_counts_55_up <- favstats(Choma5_clusts_shorter2_no_home_biting_places2_55_up$clust~Choma5_clusts_shorter2_no_home_biting_places2_55_up$partid)$n
only_home_1_16 <- nlevels(as.factor(Choma5_clusts_shorter2_1_16$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_1_16$partid))
only_home_17_34 <- nlevels(as.factor(Choma5_clusts_shorter2_17_34$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_17_34$partid))
only_home_35_54 <- nlevels(as.factor(Choma5_clusts_shorter2_35_54$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_35_54$partid))
only_home_55_up <- nlevels(as.factor(Choma5_clusts_shorter2_55_up$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_55_up$partid))
loc_counts_1_162 <- c(loc_counts_1_16,rep(0, only_home_1_16))
loc_counts_17_342 <- c(loc_counts_17_34,rep(0, only_home_17_34))
loc_counts_35_542 <- c(loc_counts_35_54,rep(0, only_home_35_54))
loc_counts_55_up2 <- c(loc_counts_55_up,rep(0, only_home_55_up))

favstats(loc_counts_1_162)
# min Q1 median Q3 max     mean       sd  n missing
#    0  0      2  3   5 1.818 1.779 11       0
favstats(loc_counts_17_342)
# min Q1 median Q3 max     mean       sd  n missing
#   0  0      1  3   6 1.933 2.219 15       0
favstats(loc_counts_35_542)
# min Q1 median Q3 max     mean       sd  n missing
#   0  1      1  3   6 2.045 1.759 22       0
favstats(loc_counts_55_up2)
# min Q1 median Q3 max     mean       sd  n missing
#   0  1      2 2.75   5 1.786 1.424 14       0

summary(as.factor(loc_counts_1_162))
# 0 1 2 3 4 5 
# 4 1 2 2 1 1 
summary(as.factor(loc_counts_17_342))
# 0 1 2 3 5 6 
# 6 2 2 2 1 2 
summary(as.factor(loc_counts_35_542))
# 0 1 2 3 4 5 6 
# 4 8 1 4 3 1 1 
summary(as.factor(loc_counts_55_up2))
# 0 1 2 3 5
# 3 3 4 3 1 
## all those with 0 have biting time at home location (not just nowhere)

locs_count2 <- as.data.frame(cbind(c(loc_counts_1_162, loc_counts_17_342, loc_counts_35_542, loc_counts_55_up2),
                                   c(rep("1-16",length(loc_counts_1_162)),
                                     rep("17-34",length(loc_counts_17_342)),
                                     rep("35-54",length(loc_counts_35_542)),
                                     rep("55+",length(loc_counts_55_up2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) + 
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-0.35,6.35)) +
  scale_x_continuous(breaks=c(seq(0,7,1)))
#
##
### time per location ####
favstats((Choma5_clusts_shorter2_no_home_biting_places2_1_16$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5636 0.884  2.079 21.8 218.4 39.47 73.18 20       0
favstats((Choma5_clusts_shorter2_no_home_biting_places2_17_34$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5383 0.9056  2.315 12.28 351.6 39.57 91.42 29       0
favstats((Choma5_clusts_shorter2_no_home_biting_places2_35_54$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5142 2.06   8.42 29.25 291.3 45.1 78.12 45       0
favstats((Choma5_clusts_shorter2_no_home_biting_places2_55_up$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.615 2.621  16.72 161.9 269.2 75.3 90.55 25       0

clusts_median_Choma_1_16 <- (favstats(Choma5_clusts_shorter2_no_home_biting_places2_1_16$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_1_16$partid)$median)*24
clusts_median_Choma_17_34 <- (favstats(Choma5_clusts_shorter2_no_home_biting_places2_17_34$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_17_34$partid)$median)*24
clusts_median_Choma_35_54 <- (favstats(Choma5_clusts_shorter2_no_home_biting_places2_35_54$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_35_54$partid)$median)*24
clusts_median_Choma_55_up <- (favstats(Choma5_clusts_shorter2_no_home_biting_places2_55_up$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_55_up$partid)$median)*24
favstats(clusts_median_Choma_1_16)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.8986 1.573  20.62 81.74 218.4 58.08 79.64 7       0
favstats(clusts_median_Choma_17_34)
# min       Q1   median  Q3      max     mean       sd  n missing
# 1.56 1.846  4.041 7.58 114.1 17.52 36.51 9       0
favstats(clusts_median_Choma_35_54)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.5928 2.781  10.93 45.94 291.3 53.53 90.23 18       0
favstats(clusts_median_Choma_55_up)
# min       Q1   median  Q3      max     mean       sd  n missing
# 2.184 5.593  92.26 143.8 175.1 78.39 73.18 11       0

clust_time_part_age <- as.data.frame(cbind(c(clusts_median_Choma_1_16, clusts_median_Choma_17_34, clusts_median_Choma_35_54, clusts_median_Choma_55_up),
                                              c(rep("1-16",length(clusts_median_Choma_1_16)),
                                                rep("17-34",length(clusts_median_Choma_17_34)),
                                                rep("35-54",length(clusts_median_Choma_35_54)),
                                                rep("55+",length(clusts_median_Choma_55_up)))))
colnames(clust_time_part_age) <- c("values","type")
clust_time_part_age$values <- as.numeric(as.character(clust_time_part_age$values))
clust_time_part_age$type <- factor(clust_time_part_age$type)

ggplot(data=clust_time_part_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=4) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-10,300)) +
  scale_x_continuous(breaks=c(seq(0,300,50)))
#
### number trips ####
Choma5_clust_trips_shorter2_no_home_biting_places3_1_16 <- Choma5_clust_trips_shorter2_no_home_biting_places3[which(Choma5_clust_trips_shorter2_no_home_biting_places3$partid %in% Choma_parts_1_16),]
Choma5_clust_trips_shorter2_no_home_biting_places3_17_34 <- Choma5_clust_trips_shorter2_no_home_biting_places3[which(Choma5_clust_trips_shorter2_no_home_biting_places3$partid %in% Choma_parts_17_34),]
Choma5_clust_trips_shorter2_no_home_biting_places3_35_54 <- Choma5_clust_trips_shorter2_no_home_biting_places3[which(Choma5_clust_trips_shorter2_no_home_biting_places3$partid %in% Choma_parts_35_54),]
Choma5_clust_trips_shorter2_no_home_biting_places3_55_up <- Choma5_clust_trips_shorter2_no_home_biting_places3[which(Choma5_clust_trips_shorter2_no_home_biting_places3$partid %in% Choma_parts_55_up),]

favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#    1  1      1  2  23 3.895 6.757 19       0
favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2  22 3.207 5.088 29       0
favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      2  4  27 4.289 6.025 45       0
favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#    1  1      1  2  22 3.207 5.088 29       0

summary(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2))
#  1  2  3  8 22 23 
# 12  3  1  1  1  1 
summary(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2))
#   1  2  3  4 12 13 22 
#  20  3  1  1  1  2  1
summary(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2))
#   1  2  3  4  5  7  8 10 12 14 15 17 22 27 
#  21 10  2  1  2  1  1  1  1  1  1  1  1  1 
summary(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2))
#  1  2  3  4  5  6  8 10 21 27 
#  9  4  4  1  1  1  1  2  1  1 

### PER PART ###
trip_count_median_1_16 <- favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_1_16$partid)$median
favstats(trip_count_median_1_16)
# min Q1 median Q3 max  mean  sd  n missing
#    1  1    1.5 3.25  22 4.714 7.724 7       0
trip_count_median_17_34 <- favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_17_34$partid)$median
favstats(trip_count_median_17_34)
# min Q1 median Q3 max  mean  sd  n missing
#  1  1      1  1   7 1.778 1.986 9       0
trip_count_median_35_54 <- favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_35_54$partid)$median
favstats(trip_count_median_35_54)
# min Q1 median Q3 max  mean  sd  n missing
#   1 1.125      2 3.875  27 4.361 6.405 18       0
trip_count_median_55_up <- favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_55_up$partid)$median
favstats(trip_count_median_55_up)
# min Q1 median Q3 max  mean  sd  n missing
#   1  2      3  4   6 3.136 1.66 11       0

trips_part2_age <- as.data.frame(cbind(c(trip_count_median_1_16, trip_count_median_17_34, trip_count_median_35_54, trip_count_median_55_up),
                                          c(rep("1-16",length(trip_count_median_1_16)),
                                            rep("17-34",length(trip_count_median_17_34)),
                                            rep("35-54",length(trip_count_median_35_54)),
                                            rep("55+",length(trip_count_median_55_up)))))
colnames(trips_part2_age) <- c("values","type")
trips_part2_age$values <- as.numeric(as.character(trips_part2_age$values))
trips_part2_age$type <- factor(trips_part2_age$type)

ggplot(data=trips_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,26.25)) + 
  scale_x_continuous(breaks=c(seq(0,30,2)))


## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,30))+ 
  scale_x_continuous(breaks=c(seq(0,30,2)))
#
##### percent biting time spent at home vs elsewhere #####
Choma5_clusts_shorter2_biting_places3b_1_16 <- Choma5_clusts_shorter2_biting_places3b[which(Choma5_clusts_shorter2_biting_places3b$partid %in% Choma_parts_1_16),]
Choma5_clusts_shorter2_biting_places3b_17_34 <- Choma5_clusts_shorter2_biting_places3b[which(Choma5_clusts_shorter2_biting_places3b$partid %in% Choma_parts_17_34),]
Choma5_clusts_shorter2_biting_places3b_35_54 <- Choma5_clusts_shorter2_biting_places3b[which(Choma5_clusts_shorter2_biting_places3b$partid %in% Choma_parts_35_54),]
Choma5_clusts_shorter2_biting_places3b_55_up <- Choma5_clusts_shorter2_biting_places3b[which(Choma5_clusts_shorter2_biting_places3b$partid %in% Choma_parts_55_up),]

favstats(Choma5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time)
favstats(Choma5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time[which(Choma5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time != 0)])
favstats(Choma5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time)
favstats(Choma5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time[which(Choma5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time != 0)])
favstats(Choma5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time)
favstats(Choma5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time[which(Choma5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time != 0)])
favstats(Choma5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time)
favstats(Choma5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time[which(Choma5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time != 0)])


perc_time_home2_age <- as.data.frame(cbind(c(Choma5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time, Choma5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time,
                                                Choma5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time, Choma5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time),
                                              c(rep("1-16",length(Choma5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time)),
                                                rep("17-34",length(Choma5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time)),
                                                rep("35-54",length(Choma5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time)),
                                                rep("55+",length(Choma5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time)))))
colnames(perc_time_home2_age) <- c("values","type")
perc_time_home2_age$values <- as.numeric(as.character(perc_time_home2_age$values))
perc_time_home2_age$type <- factor(perc_time_home2_age$type)

ggplot(data=perc_time_home2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-1,101))
##
###############
###############
######### MUTASA #########
###### General Age Analysis ####
Mutasa_part_info <- Mutasa5_clusts_shorter2[!duplicated(Mutasa5_clusts_shorter2$partid),c(1:3)]
summary(Mutasa_part_info$age_all)
# 4 - 88; mean=38.3, median=35, 1Q=26, 3Q=49.75

Mutasa_part_info$partid2 <- gsub("b","", Mutasa_part_info$partid)
Mutasa_part_info2 <- Mutasa_part_info[!duplicated(Mutasa_part_info$partid2),]
summary(as.numeric(Mutasa_part_info2$age_all))


Mutasa_part_info$age_18 <- factor(0, levels=c(0,1))
Mutasa_part_info$age_18[which(Mutasa_part_info$age_all >= 18)] <- 1
summary(as.factor(Mutasa_part_info$age_18))
#  0   1
# 27  153  
Mutasa_part_info$age_16_35_55_up <- factor(0, levels=c(0,1,2,3)) # 0 = <16; 1 = 17-34; 2 = 35-54; 3= 55+
Mutasa_part_info$age_16_35_55_up[which(Mutasa_part_info$age_all >= 17 & Mutasa_part_info$age_all < 35)] <- 1
Mutasa_part_info$age_16_35_55_up[which(Mutasa_part_info$age_all >= 35 & Mutasa_part_info$age_all < 55)] <- 2
Mutasa_part_info$age_16_35_55_up[which(Mutasa_part_info$age_all >= 55)] <- 3
summary(as.factor(Mutasa_part_info$age_16_35_55_up))
#  0   1   2   3
# 21  68  57  34
Mutasa_part_info$age_16_30_45_60_up <- factor(0, levels=c(0,1,2,3,4)) # 0 = <16; 1 = 17-29; 2 = 30-45; 3 = 46-59; 4= 60+
Mutasa_part_info$age_16_30_45_60_up[which(Mutasa_part_info$age_all >= 17 & Mutasa_part_info$age_all < 30)] <- 1
Mutasa_part_info$age_16_30_45_60_up[which(Mutasa_part_info$age_all >= 30 & Mutasa_part_info$age_all < 45)] <- 2
Mutasa_part_info$age_16_30_45_60_up[which(Mutasa_part_info$age_all >= 45 & Mutasa_part_info$age_all < 60)] <- 3
Mutasa_part_info$age_16_30_45_60_up[which(Mutasa_part_info$age_all >= 60)] <- 4
summary(as.factor(Mutasa_part_info$age_16_30_45_60_up))
#  0   1   2   3   4
# 21  39  65  33  22
ggplot(data=Mutasa_part_info, aes(x=age_all)) +
  geom_histogram(binwidth = 1, aes(fill=age_16_35_55_up)) +
  coord_cartesian(xlim=c(0,90))

ggplot(data=Mutasa_part_info, aes(x=age_all)) +
  geom_histogram(binwidth = 1, aes(fill=age_16_30_45_60_up)) +
  coord_cartesian(xlim=c(0,90))
#

###### Split by Age ######
Mutasa_parts_1_16 <- Mutasa_part_info$partid[which(Mutasa_part_info$age_16_35_55_up == 0)]
Mutasa_parts_17_34 <- Mutasa_part_info$partid[which(Mutasa_part_info$age_16_35_55_up == 1)]
Mutasa_parts_35_54 <- Mutasa_part_info$partid[which(Mutasa_part_info$age_16_35_55_up == 2)]
Mutasa_parts_55_up <- Mutasa_part_info$partid[which(Mutasa_part_info$age_16_35_55_up == 3)]
##
### number of locations (without home) ####
Mutasa5_clusts_shorter2_no_home_1_16 <- Mutasa5_clusts_shorter2_no_home[which(Mutasa5_clusts_shorter2_no_home$partid %in% Mutasa_parts_1_16),]
Mutasa5_clusts_shorter2_no_home_17_34 <- Mutasa5_clusts_shorter2_no_home[which(Mutasa5_clusts_shorter2_no_home$partid %in% Mutasa_parts_17_34),]
Mutasa5_clusts_shorter2_no_home_35_54 <- Mutasa5_clusts_shorter2_no_home[which(Mutasa5_clusts_shorter2_no_home$partid %in% Mutasa_parts_35_54),]
Mutasa5_clusts_shorter2_no_home_55_up <- Mutasa5_clusts_shorter2_no_home[which(Mutasa5_clusts_shorter2_no_home$partid %in% Mutasa_parts_55_up),]

loc_counts_1_16 <- favstats(Mutasa5_clusts_shorter2_no_home_1_16$clust~Mutasa5_clusts_shorter2_no_home_1_16$partid)$n
only_home_1_16 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_1_16)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_1_16$partid))
loc_counts_1_162 <- c(loc_counts_1_16,rep(0, only_home_1_16))
favstats(loc_counts_1_162)
# min   Q1 median     Q3  max     mean        sd  n missing
#    1  3      6  8  13 6.142857 3.863751 21       0
loc_counts_17_34 <- favstats(Mutasa5_clusts_shorter2_no_home_17_34$clust~Mutasa5_clusts_shorter2_no_home_17_34$partid)$n
only_home_17_34 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_17_34)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_17_34$partid))
loc_counts_17_342 <- c(loc_counts_17_34,rep(0, only_home_17_34))
favstats(loc_counts_17_342)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0  4      6 11  25 7.971014 6.053611 69       0
loc_counts_35_54 <- favstats(Mutasa5_clusts_shorter2_no_home_35_54$clust~Mutasa5_clusts_shorter2_no_home_35_54$partid)$n
only_home_35_54 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_35_54)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_35_54$partid))
loc_counts_35_542 <- c(loc_counts_35_54,rep(0, only_home_35_54))
favstats(loc_counts_35_542)
# min   Q1 median     Q3  max     mean        sd  n missing
#     0  4    7.5 12  21 8.396552 5.681425 58       0
loc_counts_55_up <- favstats(Mutasa5_clusts_shorter2_no_home_55_up$clust~Mutasa5_clusts_shorter2_no_home_55_up$partid)$n
only_home_55_up <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_55_up)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_55_up$partid))
loc_counts_55_up2 <- c(loc_counts_55_up,rep(0, only_home_55_up))
favstats(loc_counts_55_up2)
# min   Q1 median     Q3  max     mean        sd  n missing
#  1  4    6.5 11.75  24 8.176471 5.926165 34       0

agecols <- brewer.pal(5, "YlGn")

locs_count2 <- as.data.frame(cbind(c(loc_counts_1_162, loc_counts_17_342, loc_counts_35_542, loc_counts_55_up2),
                                   c(rep("1-16",length(loc_counts_1_162)),
                                     rep("17-34",length(loc_counts_17_342)),
                                     rep("35-54",length(loc_counts_35_542)),
                                     rep("55+",length(loc_counts_55_up2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

shapiro.test(locs_count2$values) ## not normal
kruskal.test(locs_count2$values, locs_count2$type) # p = 0.54
pairwise.wilcox.test(locs_count2$values, locs_count2$type, p.adjust.method = "holm")

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-2,27)) +
  scale_x_continuous(breaks=c(seq(0,27,5)))
#
### distance of locations (without home) ####
### TOTAL ###
favstats(Mutasa5_clusts_shorter2_no_home_1_16$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
# 0.0350883 0.6723648 1.348071 2.667145 50.21874 2.93685 6.447607 129       0
favstats(Mutasa5_clusts_shorter2_no_home_17_34$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 2.745604e-05 0.5252042 1.438887 4.606841 226.2527 17.3314 49.50135 550       0
favstats(Mutasa5_clusts_shorter2_no_home_35_54$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#  0.006409266 0.5262203 1.411451 8.004257 481.4493 27.23318 74.67346 487       0
favstats(Mutasa5_clusts_shorter2_no_home_55_up$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 0.03568206 0.3530347 0.9854358 2.421828 213.6842 4.541411 20.65232 278       0

### PER PART ###
medians_km_1_16 <- (favstats(Mutasa5_clusts_shorter2_no_home_1_16$hhdist_m_haversine~Mutasa5_clusts_shorter2_no_home_1_16$partid)$median)/1000
medians_km_17_34 <- (favstats(Mutasa5_clusts_shorter2_no_home_17_34$hhdist_m_haversine~Mutasa5_clusts_shorter2_no_home_17_34$partid)$median)/1000
medians_km_35_54 <- (favstats(Mutasa5_clusts_shorter2_no_home_35_54$hhdist_m_haversine~Mutasa5_clusts_shorter2_no_home_35_54$partid)$median)/1000
medians_km_55_up <- (favstats(Mutasa5_clusts_shorter2_no_home_55_up$hhdist_m_haversine~Mutasa5_clusts_shorter2_no_home_55_up$partid)$median)/1000

favstats(medians_km_1_16) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.08078032 0.8409646 1.234995 2.033209 5.392751 1.665432 1.354737 21       0
favstats(medians_km_17_34) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.1094686 0.627445 1.297847 3.166249 217.5585 10.50969 36.66726 67       0
favstats(medians_km_35_54) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.1665696 0.8255315 1.452268 8.617983 299.3132 19.72807 52.69356 54       0
favstats(medians_km_55_up) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.1013454 0.5272192 0.8293977 1.381213 6.551693 1.166241 1.185106 34       0

medians_km_17_34_short <- medians_km_17_34[-which(medians_km_17_34 > 100)]
medians_km_35_54_short <- medians_km_35_54[-which(medians_km_35_54 > 100)]
favstats(medians_km_17_34_short) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
# 0.1094686 0.6262455 1.296415 3.102046 42.06424 4.300304 8.591922 65       0
favstats(medians_km_35_54_short) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
# 0.1665696 0.7659092 1.278861 2.692459 57.42812 6.76014 13.1404 50       0

dists_km2_age <- as.data.frame(cbind(c(medians_km_1_16, medians_km_17_34_short, medians_km_35_54_short, medians_km_55_up),
                                     c(rep("1-16",length(medians_km_1_16)),
                                       rep("17-34",length(medians_km_17_34_short)),
                                       rep("35-54",length(medians_km_35_54_short)),
                                       rep("55+",length(medians_km_55_up)))))
colnames(dists_km2_age) <- c("values","type")
dists_km2_age$values <- as.numeric(as.character(dists_km2_age$values))
dists_km2_age$type <- factor(dists_km2_age$type)

shapiro.test(dists_km2_age$values) ## not normal
kruskal.test(dists_km2_age$values, dists_km2_age$type) # p = 0.04
pairwise.wilcox.test(dists_km2_age$values, dists_km2_age$type, p.adjust.method = "holm")

ggplot(data=dists_km2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=12) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-3,65))
###
### number of trips (without home) ####
### TOTAL ###
Mutasa5_clust_trips_shorter2_no_home3_1_16 <- Mutasa5_clust_trips_shorter2_no_home3[which(Mutasa5_clust_trips_shorter2_no_home3$partid %in% Mutasa_parts_1_16),]
Mutasa5_clust_trips_shorter2_no_home3_17_34 <- Mutasa5_clust_trips_shorter2_no_home3[which(Mutasa5_clust_trips_shorter2_no_home3$partid %in% Mutasa_parts_17_34),]
Mutasa5_clust_trips_shorter2_no_home3_35_54 <- Mutasa5_clust_trips_shorter2_no_home3[which(Mutasa5_clust_trips_shorter2_no_home3$partid %in% Mutasa_parts_35_54),]
Mutasa5_clust_trips_shorter2_no_home3_55_up <- Mutasa5_clust_trips_shorter2_no_home3[which(Mutasa5_clust_trips_shorter2_no_home3$partid %in% Mutasa_parts_55_up),]

favstats(Mutasa5_clust_trips_shorter2_no_home3_1_16$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      2  4  29 2.930233 3.750388 129       0
favstats(Mutasa5_clust_trips_shorter2_no_home3_17_34$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      1  2  39 2.501818 3.708283 550       0
favstats(Mutasa5_clust_trips_shorter2_no_home3_35_54$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      1  3  34 2.533881 3.165916 487       0
favstats(Mutasa5_clust_trips_shorter2_no_home3_55_up$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
# 1  1      2  3  43 2.658273 4.056699 278       0

### PER PART ###
trips_part_median_1_16 <- favstats(Mutasa5_clust_trips_shorter2_no_home3_1_16$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_1_16$partid)$median
trips_part_median_17_34 <- favstats(Mutasa5_clust_trips_shorter2_no_home3_17_34$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_17_34$partid)$median
trips_part_median_35_54 <- favstats(Mutasa5_clust_trips_shorter2_no_home3_35_54$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_35_54$partid)$median
trips_part_median_55_up <- favstats(Mutasa5_clust_trips_shorter2_no_home3_55_up$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_55_up$partid)$median

favstats(trips_part_median_1_16)
#   min  Q1   median   Q3      max       mean       sd  n missing
#     1  1      2  2  13 2.404762 2.567192 21       0
favstats(trips_part_median_17_34)
#   min Q1 median  Q3   max      mean       sd  n missing
#    1  1    1.5  2  13 1.723881 1.54793 67       0
favstats(trips_part_median_35_54)
#   min  Q1   median   Q3      max       mean       sd  n missing
#     1  1    1.5  2   7 1.666667 0.9615239 54       0
favstats(trips_part_median_55_up)
#   min Q1 median  Q3   max      mean       sd  n missing
#    1  1      2 2.375 11.5 2.264706 1.985689 34       0

trips_part2_age <- as.data.frame(cbind(c(trips_part_median_1_16, trips_part_median_17_34, trips_part_median_35_54, trips_part_median_55_up),
                                       c(rep("1-16",length(trips_part_median_1_16)),
                                         rep("17-34",length(trips_part_median_17_34)),
                                         rep("35-54",length(trips_part_median_35_54)),
                                         rep("55+",length(trips_part_median_55_up)))))
colnames(trips_part2_age) <- c("values","type")
trips_part2_age$values <- as.numeric(as.character(trips_part2_age$values))
trips_part2_age$type <- factor(trips_part2_age$type)
kruskal.test(trips_part2_age$values, trips_part2_age$type) # p = 0.04
pairwise.wilcox.test(trips_part2_age$values, trips_part2_age$type, p.adjust.method = "holm")

ggplot(data=trips_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,13)) +
  scale_x_continuous(breaks=c(0:13))
#
#
### time per trip ####
### TOTAL ###
Mutasa5_clust_trips_shorter2_no_home_1_16 <- Mutasa5_clust_trips_shorter2_no_home[which(Mutasa5_clust_trips_shorter2_no_home$partid %in% Mutasa_parts_1_16),]
Mutasa5_clust_trips_shorter2_no_home_17_34 <- Mutasa5_clust_trips_shorter2_no_home[which(Mutasa5_clust_trips_shorter2_no_home$partid %in% Mutasa_parts_17_34),]
Mutasa5_clust_trips_shorter2_no_home_35_54 <- Mutasa5_clust_trips_shorter2_no_home[which(Mutasa5_clust_trips_shorter2_no_home$partid %in% Mutasa_parts_35_54),]
Mutasa5_clust_trips_shorter2_no_home_55_up <- Mutasa5_clust_trips_shorter2_no_home[which(Mutasa5_clust_trips_shorter2_no_home$partid %in% Mutasa_parts_55_up),]

favstats(Mutasa5_clust_trips_shorter2_no_home_1_16$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.48 1.531458 2.747222 6.262222 162.2475 5.697354 13.8073 378       0
favstats(Mutasa5_clust_trips_shorter2_no_home_17_34$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.4413889 1.271806 2.400278 4.761667 1218.114 9.217612 53.2231 1376       0
favstats(Mutasa5_clust_trips_shorter2_no_home_35_54$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.48 1.436597 2.613333 5.674444 331.2 6.40675 19.44737 1234       0
favstats(Mutasa5_clust_trips_shorter2_no_home_55_up$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.5172222 1.494722 2.742222 5.210278 58.45806 4.017645 4.672734 739       0

### PER PART ###
trip_time_part_median_1_16 <- favstats(Mutasa5_clust_trips_shorter2_no_home_1_16$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_1_16$partid)$median
trip_time_part_median_17_34 <- favstats(Mutasa5_clust_trips_shorter2_no_home_17_34$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_17_34$partid)$median
trip_time_part_median_35_54 <- favstats(Mutasa5_clust_trips_shorter2_no_home_35_54$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_35_54$partid)$median
trip_time_part_median_55_up <- favstats(Mutasa5_clust_trips_shorter2_no_home_55_up$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_55_up$partid)$median

favstats(trip_time_part_median_1_16)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 1.506389 2.009028 3.801944 6.245833 9.897222 4.354577 2.571019 21       0
favstats(trip_time_part_median_17_34)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.9602778 2.031389 2.728889 3.49875 1218.114 21.21383 148.451 67       0
favstats(trip_time_part_median_35_54)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 1.440139 2.229861 2.859444 4.091667 33.96292 4.181847 4.724037 54       0
favstats(trip_time_part_median_55_up)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.192639 1.972847 2.641528 3.341736 7.917222 2.97011 1.456799 34       0

## REMOVE OUTLIER ##
trip_time_part_median_17_34_short <- trip_time_part_median_17_34[-which(trip_time_part_median_17_34 > 600)]
favstats(trip_time_part_median_17_34_short)
# min       Q1   median       Q3      max     mean       sd   n missing
# 0.9602778 2.030556 2.719583 3.450486 12.22667 3.078984 1.810434 66       0
trip_time_part_median_35_54_short <- trip_time_part_median_35_54[-which(trip_time_part_median_35_54 > 30)]
favstats(trip_time_part_median_35_54_short)
# min       Q1   median       Q3      max     mean       sd   n missing
# 1.440139 2.227222 2.839167 3.995417 12.13694 3.61994 2.316868 53       0



trip_times_part2_age <- as.data.frame(cbind(c(trip_time_part_median_1_16, trip_time_part_median_17_34_short, trip_time_part_median_35_54_short, trip_time_part_median_55_up),
                                            c(rep("1-16",length(trip_time_part_median_1_16)),
                                              rep("17-34",length(trip_time_part_median_17_34_short)),
                                              rep("35-54",length(trip_time_part_median_35_54_short)),
                                              rep("55+",length(trip_time_part_median_55_up)))))
colnames(trip_times_part2_age) <- c("values","type")
trip_times_part2_age$values <- as.numeric(as.character(trip_times_part2_age$values))
trip_times_part2_age$type <- factor(trip_times_part2_age$type)
kruskal.test(trip_times_part2_age$values, trip_times_part2_age$type) # p = 0.04
pairwise.wilcox.test(trip_times_part2_age$values, trip_times_part2_age$type, p.adjust.method = "holm")

ggplot(data=trip_times_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-0.5,13)) +
  scale_x_continuous(breaks=c(seq(0,14,2)))
#
#

### Home Loc ####
Mutasa5_clusts_shorter2_home_1_16 <- Mutasa5_clusts_shorter2_home[which(Mutasa5_clusts_shorter2_home$partid %in% Mutasa_parts_1_16),] 
no_home_1_16 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_1_16)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_home_1_16$partid))# 1
## none with no home
Mutasa5_clusts_shorter2_home_17_34 <- Mutasa5_clusts_shorter2_home[which(Mutasa5_clusts_shorter2_home$partid %in% Mutasa_parts_17_34),] 
no_home_17_34 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_17_34)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_home_17_34$partid))# 1
## none with no home
Mutasa5_clusts_shorter2_home_35_54 <- Mutasa5_clusts_shorter2_home[which(Mutasa5_clusts_shorter2_home$partid %in% Mutasa_parts_35_54),] 
no_home_35_54 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_35_54)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_home_35_54$partid))# 1
## none with no home
Mutasa5_clusts_shorter2_home_55_up <- Mutasa5_clusts_shorter2_home[which(Mutasa5_clusts_shorter2_home$partid %in% Mutasa_parts_55_up),] 
no_home_55_up <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_55_up)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_home_55_up$partid))# 1
## none with no home

### number of 'trips' to home loc ####
### TOTAL ###
Mutasa5_clust_trips_shorter2_home3_1_16 <- Mutasa5_clust_trips_shorter2_home3[which(Mutasa5_clust_trips_shorter2_home3$partid %in% Mutasa_parts_1_16),]
Mutasa5_clust_trips_shorter2_home3_17_34 <- Mutasa5_clust_trips_shorter2_home3[which(Mutasa5_clust_trips_shorter2_home3$partid %in% Mutasa_parts_17_34),]
Mutasa5_clust_trips_shorter2_home3_35_54 <- Mutasa5_clust_trips_shorter2_home3[which(Mutasa5_clust_trips_shorter2_home3$partid %in% Mutasa_parts_35_54),]
Mutasa5_clust_trips_shorter2_home3_55_up <- Mutasa5_clust_trips_shorter2_home3[which(Mutasa5_clust_trips_shorter2_home3$partid %in% Mutasa_parts_55_up),]

favstats(Mutasa5_clust_trips_shorter2_home3_1_16$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#      2  8     11 18  31 13.28571 8.355494 21       0
favstats(Mutasa5_clust_trips_shorter2_home3_17_34$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#     1  4      7 12.5  31 9.676471 7.769902 68       0
favstats(Mutasa5_clust_trips_shorter2_home3_35_54$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#   1  4     12 19  40 12.80702 9.9093 57       0
favstats(Mutasa5_clust_trips_shorter2_home3_55_up$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#    3  7     25 23  42 17.73529 11.13365 34       0
#
### total time at home loc ####
Mutasa5_clust_trips_shorter2_home4_1_16 <- Mutasa5_clust_trips_shorter2_home4[which(Mutasa5_clust_trips_shorter2_home4$partid %in% Mutasa_parts_1_16),]
favstats(Mutasa5_clust_trips_shorter2_home4_1_16$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
#  6.948924 21.85723 25.013 33.86895 60.75051 27.29022 12.03058 21       0
Mutasa5_clust_trips_shorter2_home4_17_34 <- Mutasa5_clust_trips_shorter2_home4[which(Mutasa5_clust_trips_shorter2_home4$partid %in% Mutasa_parts_17_34),]
favstats(Mutasa5_clust_trips_shorter2_home4_17_34$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
#  0.925081 15.92943 24.81459 28.77078 54.91443 24.06588 11.45985 68       0
Mutasa5_clust_trips_shorter2_home4_35_54 <- Mutasa5_clust_trips_shorter2_home4[which(Mutasa5_clust_trips_shorter2_home4$partid %in% Mutasa_parts_35_54),]
favstats(Mutasa5_clust_trips_shorter2_home4_35_54$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
# 0.02479167 17.84415 21.37396 32.53422 65.96275 24.82486 11.94777 57       0
Mutasa5_clust_trips_shorter2_home4_55_up <- Mutasa5_clust_trips_shorter2_home4[which(Mutasa5_clust_trips_shorter2_home4$partid %in% Mutasa_parts_55_up),]
favstats(Mutasa5_clust_trips_shorter2_home4_55_up$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
#  8.41285 15.65769 23.09618 27.23998 59.14673 23.67459 11.11103 34       0
#
### percent time at home ####
favstats(Mutasa5_clusts_shorter2_home_1_16$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
#  50.64093 79.88855 90.28389 93.49001 99.38704 85.58539 12.22839 21       0
favstats(Mutasa5_clusts_shorter2_home_17_34$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#  14.3927 65.87306 87.06052 97.31496 100 78.14987 23.48485 69       0
favstats(Mutasa5_clusts_shorter2_home_35_54$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
# 31.3768 74.96783 88.17172 95.05399 100 82.74607 16.5383 58       0
favstats(Mutasa5_clusts_shorter2_home_55_up$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#  55.14049 73.76317 88.22443 97.07727 99.79577 85.19243 12.61589 34       0

perc_home2_age <- as.data.frame(cbind(c(Mutasa5_clusts_shorter2_home_1_16$perc_time_home, Mutasa5_clusts_shorter2_home_17_34$perc_time_home, Mutasa5_clusts_shorter2_home_35_54$perc_time_home, Mutasa5_clusts_shorter2_home_55_up$perc_time_home),
                                      c(rep("1-16",length(Mutasa5_clusts_shorter2_home_1_16$perc_time_home)),
                                        rep("17-34",length(Mutasa5_clusts_shorter2_home_17_34$perc_time_home)),
                                        rep("35-54",length(Mutasa5_clusts_shorter2_home_35_54$perc_time_home)),
                                        rep("55+",length(Mutasa5_clusts_shorter2_home_55_up$perc_time_home)))))
colnames(perc_home2_age) <- c("values","type")
perc_home2_age$values <- as.numeric(as.character(perc_home2_age$values))
perc_home2_age$type <- factor(perc_home2_age$type)

kruskal.test(perc_home2_age$values, perc_home2_age$type) # p = 0.04
pairwise.wilcox.test(perc_home2_age$values, perc_home2_age$type, p.adjust.method = "holm")

ggplot(data=perc_home2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,101))
##
### Average Direction of Travel #####
temp <- Mutasa5_clusts_shorter2_no_home[!duplicated(Mutasa5_clusts_shorter2_no_home$partid),]

age_group <- Mutasa5_clusts_shorter2_no_home$age_all[!duplicated(Mutasa5_clusts_shorter2_no_home$partid)]

temp$age_group <- factor(0, levels=c(3,2,1,0)) # 0 = <16; 1 = 17-34; 2 = 35-54; 3= 55+
temp$age_group[which(age_group >= 17 & age_group < 35)] <- 1
temp$age_group[which(age_group >= 35 & age_group < 55)] <- 2
temp$age_group[which(age_group >= 55)] <- 3
summary(as.factor(temp$age_group))

ggplot(data=temp, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=age_group), position=position_dodge2(preserve = "single"), width=4) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Age Group", values=agecols[c(5,4,2,1)], labels=c("55+", "35-54", "17-34", "0-16"))
#

favstats(temp$avg_dir_from_home_bearing~temp$age_group)
# temp$age_group       min       Q1   median       Q3      max     mean        sd  n missing
# 1              3  8.648831 127.8731 225.3825 283.1798 359.1797 250.7729  98.45838 34       0
# 2              2  6.865241 126.0778 227.8045 283.0978 342.0978 255.9048  98.49494 54       0
# 3              1 16.387774 118.5747 255.8879 273.8240 359.7089 194.9276  97.21048 67       0
# 4              0  2.636555 180.5733 242.2535 318.8992 345.2589 222.7989 117.69572 21       0
favstats(temp$avg_dir_from_home_mag~temp$age_group)
# temp$age_group        min        Q1    median        Q3        max      mean        sd  n
# 1              3 0.05659092 0.3493925 0.7625568  2.380696  29.086386  2.474395  5.087821 34
# 2              2 0.07662500 0.7162572 1.8375864 19.427715 181.957304 22.975150 44.393275 54
# 3              1 0.08965062 0.7159847 1.6909944  8.831621 210.008021 14.789211 37.548565 67
# 4              0 0.03380367 0.7047252 0.8848382  2.254216   9.488883  2.044065  2.609356 21

temp2 <- temp[which(temp$avg_dir_from_home_mag < 70),] ## 164 of 176 parts

favstats(temp2$avg_dir_from_home_bearing~temp2$age_group)
# temp2$age_group       min       Q1   median       Q3      max     mean        sd  n
# 1               3  8.648831 127.8731 225.3825 283.1798 359.1797 250.7729  98.45838 34
# 2               2  6.865241 112.4042 225.6456 268.1182 342.0978 252.3748  99.31925 47
# 3               1 16.387774 107.6938 252.0140 264.8290 359.7089 187.7596  97.57417 62
# 4               0  2.636555 180.5733 242.2535 318.8992 345.2589 222.7989 117.69572 21
favstats(temp2$avg_dir_from_home_mag~temp2$age_group)
# temp2$age_group        min        Q1    median       Q3       max     mean        sd  n
# 1               3 0.05659092 0.3493925 0.7625568 2.380696 29.086386 2.474395  5.087821 34
# 2               2 0.07662500 0.5827052 1.5216457 7.296525 68.506458 7.794781 15.019433 47
# 3               1 0.08965062 0.6348413 1.4075333 4.387003 58.354660 5.674589 10.845536 62
# 4               0 0.03380367 0.7047252 0.8848382 2.254216  9.488883 2.044065  2.609356 21

kruskal.test(temp2$avg_dir_from_home_bearing, temp2$age_group) # p = 0.04
pairwise.wilcox.test(temp2$avg_dir_from_home_bearing, temp2$age_group, p.adjust.method = "holm")


ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=age_group), position=position_dodge2(preserve = "single"), width=4) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Age Group", values=agecols[c(5,4,2,1)], labels=c("55+", "35-54", "17-34", "0-16"))



temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)] <- -(360-temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)])
ggplot(data=temp, aes(x=avg_dir_from_home_values, y=..count.., fill=age_group)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~age_group, labeller = labeller(age_group = c("0"="0-16", "1"="17-34","2"="35-54", "3"="55+")))+
  scale_fill_manual(values=agecols[c(1,2,4,5)]) +
  scale_y_continuous(breaks=c(seq(0,8,2)))

#
### time spent at clusters vs elsewhere (travel, less important places, etc) ####
# percent time in clusters
Mutasa5_clusts_shorter2_1_16 <- Mutasa5_clusts_shorter2[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_1_16),]
Mutasa5_clusts_shorter2_17_34 <- Mutasa5_clusts_shorter2[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_17_34),]
Mutasa5_clusts_shorter2_35_54 <- Mutasa5_clusts_shorter2[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_35_54),]
Mutasa5_clusts_shorter2_55_up <- Mutasa5_clusts_shorter2[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_55_up),]

Mutasa5_clusts_shorter2_1_16$perc_time_in_all_clusts <- Mutasa5_clusts_shorter2_1_16$part_trips_time_new/Mutasa5_clusts_shorter2_1_16$total_time
Mutasa5_clusts_shorter2_17_34$perc_time_in_all_clusts <- Mutasa5_clusts_shorter2_17_34$part_trips_time_new/Mutasa5_clusts_shorter2_17_34$total_time
Mutasa5_clusts_shorter2_35_54$perc_time_in_all_clusts <- Mutasa5_clusts_shorter2_35_54$part_trips_time_new/Mutasa5_clusts_shorter2_35_54$total_time
Mutasa5_clusts_shorter2_55_up$perc_time_in_all_clusts <- Mutasa5_clusts_shorter2_55_up$part_trips_time_new/Mutasa5_clusts_shorter2_55_up$total_time

summary(Mutasa5_clusts_shorter2_1_16$perc_time_in_all_clusts[!duplicated(Mutasa5_clusts_shorter2_1_16$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8451  0.8610  0.9074  0.9118  0.9541  0.9932 
summary(Mutasa5_clusts_shorter2_17_34$perc_time_in_all_clusts[!duplicated(Mutasa5_clusts_shorter2_17_34$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8663  0.9148  0.9509  0.9451  0.9895  0.9955 
summary(Mutasa5_clusts_shorter2_35_54$perc_time_in_all_clusts[!duplicated(Mutasa5_clusts_shorter2_35_54$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6629  0.8749  0.9248  0.9012  0.9642  0.9984 
summary(Mutasa5_clusts_shorter2_55_up$perc_time_in_all_clusts[!duplicated(Mutasa5_clusts_shorter2_55_up$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7781  0.8729  0.9317  0.9195  0.9782  0.9950 

# time outside clusters
summary(Mutasa5_clusts_shorter2_1_16$total_time[!duplicated(Mutasa5_clusts_shorter2_1_16$partid)] - Mutasa5_clusts_shorter2_1_16$part_trips_time_new[!duplicated(Mutasa5_clusts_shorter2_1_16$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2509  1.2039  2.9537  2.6542  4.3409  4.8027 
summary(Mutasa5_clusts_shorter2_17_34$total_time[!duplicated(Mutasa5_clusts_shorter2_17_34$partid)] - Mutasa5_clusts_shorter2_17_34$part_trips_time_new[!duplicated(Mutasa5_clusts_shorter2_17_34$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1395  0.4087  1.4091  1.8793  3.0862  4.7086 
summary(Mutasa5_clusts_shorter2_35_54$total_time[!duplicated(Mutasa5_clusts_shorter2_35_54$partid)] - Mutasa5_clusts_shorter2_35_54$part_trips_time_new[!duplicated(Mutasa5_clusts_shorter2_35_54$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2275  1.1145  2.0146  3.3552  4.5656  9.8070 
summary(Mutasa5_clusts_shorter2_55_up$total_time[!duplicated(Mutasa5_clusts_shorter2_55_up$partid)] - Mutasa5_clusts_shorter2_55_up$part_trips_time_new[!duplicated(Mutasa5_clusts_shorter2_55_up$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2383  0.7049  2.3483  2.7655  4.5573  7.0787 
#
### Time at locations ####
favstats(Mutasa5_clusts_shorter2_no_home_1_16$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.5933333 2.170278 4.037778 10.08028 472.5967 16.69457 48.38371 129       0
favstats(Mutasa5_clusts_shorter2_no_home_17_34$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.4802778 2.340833 4.0075 8.366111 1537.291 23.06079 105.896 550       0
favstats(Mutasa5_clusts_shorter2_no_home_35_54$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.5072222 2.170833 4.234167 9.221944 570.4803 16.23394 50.22989 487       0
favstats(Mutasa5_clusts_shorter2_no_home_55_up$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.5394444 2.380069   4.62 9.450069 249.3431 10.68 24.26822 278       0

### PER PART ###
clust_time_part_median_1_16 <- (favstats(Mutasa5_clusts_shorter2_no_home_1_16$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_1_16$partid)$median)*24
clust_time_part_median_17_34 <- (favstats(Mutasa5_clusts_shorter2_no_home_17_34$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_17_34$partid)$median)*24
clust_time_part_median_35_54 <- (favstats(Mutasa5_clusts_shorter2_no_home_35_54$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_35_54$partid)$median)*24
clust_time_part_median_55_up <- (favstats(Mutasa5_clusts_shorter2_no_home_55_up$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_55_up$partid)$median)*24

favstats(clust_time_part_median_1_16)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 2.249444 3.36 4.204306 6.695833 110.4386 11.69185 24.05832 21       0
favstats(clust_time_part_median_17_34)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.433611 3.159306   3.98 5.811944 1218.114 26.12027 149.7581 67       0
favstats(clust_time_part_median_35_54)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 1.666944 3.624444 4.376944 5.459931 98.0725 7.7734 15.02491 54       0
favstats(clust_time_part_median_55_up)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.656389 3.288681 4.227083 6.904861 87.48472 8.135204 14.75239 34       0
clust_time_part_median_17_342 <- clust_time_part_median_17_34[which(clust_time_part_median_17_34 < 150)]
favstats(clust_time_part_median_17_342)
# min       Q1   median       Q3      max     mean      sd  n missing
# 1.433611 3.156111 3.960278 5.620556 48.80806 5.177797 5.87062 65       0

clust_time_part_age <- as.data.frame(cbind(c(clust_time_part_median_1_16, clust_time_part_median_17_342,clust_time_part_median_35_54,clust_time_part_median_55_up),
                                           c(rep("1-16",length(clust_time_part_median_1_16)),
                                             rep("17-34",length(clust_time_part_median_17_342)),
                                             rep("35-54",length(clust_time_part_median_35_54)),
                                             rep("55+",length(clust_time_part_median_55_up)))))
colnames(clust_time_part_age) <- c("values","type")
clust_time_part_age$values <- as.numeric(as.character(clust_time_part_age$values))
clust_time_part_age$type <- factor(clust_time_part_age$type)
kruskal.test(clust_time_part_age$values, clust_time_part_age$type) # p = 0.04
pairwise.wilcox.test(clust_time_part_age$values, clust_time_part_age$type, p.adjust.method = "holm")

ggplot(data=clust_time_part_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=7) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,110)) +
  scale_x_continuous(breaks=c(seq(0,120,10)))
#
#
#
####### Biting time metrics #####
### number locations outside home ####
Mutasa5_clusts_shorter2_no_home_biting_places2_1_16 <- Mutasa5_clusts_shorter2_no_home_biting_places2[which(Mutasa5_clusts_shorter2_no_home_biting_places2$partid %in% Mutasa_parts_1_16),]
Mutasa5_clusts_shorter2_no_home_biting_places2_17_34 <- Mutasa5_clusts_shorter2_no_home_biting_places2[which(Mutasa5_clusts_shorter2_no_home_biting_places2$partid %in% Mutasa_parts_17_34),]
Mutasa5_clusts_shorter2_no_home_biting_places2_35_54 <- Mutasa5_clusts_shorter2_no_home_biting_places2[which(Mutasa5_clusts_shorter2_no_home_biting_places2$partid %in% Mutasa_parts_35_54),]
Mutasa5_clusts_shorter2_no_home_biting_places2_55_up <- Mutasa5_clusts_shorter2_no_home_biting_places2[which(Mutasa5_clusts_shorter2_no_home_biting_places2$partid %in% Mutasa_parts_55_up),]

loc_counts_1_16 <- favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_1_16$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_1_16$partid)$n
loc_counts_17_34 <- favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_17_34$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_17_34$partid)$n
loc_counts_35_54 <- favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_35_54$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_35_54$partid)$n
loc_counts_55_up <- favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_55_up$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_55_up$partid)$n
only_home_1_16 <- nlevels(as.factor(Mutasa5_clusts_shorter2_1_16$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_1_16$partid))
only_home_17_34 <- nlevels(as.factor(Mutasa5_clusts_shorter2_17_34$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_17_34$partid))
only_home_35_54 <- nlevels(as.factor(Mutasa5_clusts_shorter2_35_54$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_35_54$partid))
only_home_55_up <- nlevels(as.factor(Mutasa5_clusts_shorter2_55_up$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_55_up$partid))
loc_counts_1_162 <- c(loc_counts_1_16,rep(0, only_home_1_16))
loc_counts_17_342 <- c(loc_counts_17_34,rep(0, only_home_17_34))
loc_counts_35_542 <- c(loc_counts_35_54,rep(0, only_home_35_54))
loc_counts_55_up2 <- c(loc_counts_55_up,rep(0, only_home_55_up))

favstats(loc_counts_1_162)
# min Q1 median Q3 max     mean       sd  n missing
#    0  0      0  1   4 0.7619 1.136 21       0
favstats(loc_counts_17_342)
# min Q1 median Q3 max     mean       sd  n missing
#    0  0      1  2   8 1.681 2.018 69       0
favstats(loc_counts_35_542)
# min Q1 median Q3 max     mean       sd  n missing
#   0  0      1  2   8 1.483 1.828 58       0
favstats(loc_counts_55_up2)
# min Q1 median Q3 max     mean       sd  n missing
#   0  0      1  2   5 1.029 1.243 34       0

summary(as.factor(loc_counts_1_162))
# 0  1  2  3  4 
# 12  5  2  1  1 
summary(as.factor(loc_counts_17_342))
# 0  1  2  3  4  6  8 
# 26 12 16  4  6  2  3 
summary(as.factor(loc_counts_35_542))
# 0  1  2  3  4  6  7  8 
# 20 19  7  6  2  2  1  1 
summary(as.factor(loc_counts_55_up2))
# 0  1  2  3  5 
# 16  7  7  3  1 
## all those with 0 have biting time at home location (not just nowhere)

locs_count2 <- as.data.frame(cbind(c(loc_counts_1_162, loc_counts_17_342, loc_counts_35_542, loc_counts_55_up2),
                                   c(rep("1-16",length(loc_counts_1_162)),
                                     rep("17-34",length(loc_counts_17_342)),
                                     rep("35-54",length(loc_counts_35_542)),
                                     rep("55+",length(loc_counts_55_up2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) + 
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-0.35,8.35)) +
  scale_x_continuous(breaks=c(seq(0,8,1)))
#
##
### time per location ####
favstats((Mutasa5_clusts_shorter2_no_home_biting_places2_1_16$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 1.278 1.552   7.44 14.05 84.88 15.7 24.56 16       0
favstats((Mutasa5_clusts_shorter2_no_home_biting_places2_17_34$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5072 1.51  6.717 17.69 504.1 25.23 60.56 116       0
favstats((Mutasa5_clusts_shorter2_no_home_biting_places2_35_54$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5274 2.21  8.197 17.21 196.4 22.21 40.08 86       0
favstats((Mutasa5_clusts_shorter2_no_home_biting_places2_55_up$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5961 3.956  8.425 9.342 33.15 8.235 6.373 35       0

clusts_median_Mutasa_1_16 <- (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_1_16$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_1_16$partid)$median)*24
clusts_median_Mutasa_17_34 <- (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_17_34$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_17_34$partid)$median)*24
clusts_median_Mutasa_35_54 <- (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_35_54$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_35_54$partid)$median)*24
clusts_median_Mutasa_55_up <- (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_55_up$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_55_up$partid)$median)*24
favstats(clusts_median_Mutasa_1_16)
# min       Q1   median  Q3      max     mean       sd  n missing
#  1.278 1.567   6.04 17.42 84.88 18.37 27.61 9       0
favstats(clusts_median_Mutasa_17_34)
# min       Q1   median  Q3      max     mean       sd  n missing
# 1.08 4.144  8.709 36.12 504.1 36.87 84.45 43       0
favstats(clusts_median_Mutasa_35_54)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.5831 3.88   8.49 19.09 172.1 20.24 35.44 38       0
favstats(clusts_median_Mutasa_55_up)
# min       Q1   median  Q3      max     mean       sd  n missing
# 1.374 5.089   8.14 8.839 17.4 7.198 3.747 18       0

clust_time_part_age <- as.data.frame(cbind(c(clusts_median_Mutasa_1_16, clusts_median_Mutasa_17_34, clusts_median_Mutasa_35_54, clusts_median_Mutasa_55_up),
                                           c(rep("1-16",length(clusts_median_Mutasa_1_16)),
                                             rep("17-34",length(clusts_median_Mutasa_17_34)),
                                             rep("35-54",length(clusts_median_Mutasa_35_54)),
                                             rep("55+",length(clusts_median_Mutasa_55_up)))))
colnames(clust_time_part_age) <- c("values","type")
clust_time_part_age$values <- as.numeric(as.character(clust_time_part_age$values))
clust_time_part_age$type <- factor(clust_time_part_age$type)

ggplot(data=clust_time_part_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=4) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-10,500)) +
  scale_x_continuous(breaks=c(seq(0,500,50)))
#

clust_time_part_age2 <- clust_time_part_age[which(clust_time_part_age$values < 500),]
ggplot(data=clust_time_part_age2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=4) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-10,250)) +
  scale_x_continuous(breaks=c(seq(0,250,50)))
#
### number trips ####
Mutasa5_clust_trips_shorter2_no_home_biting_places3_1_16 <- Mutasa5_clust_trips_shorter2_no_home_biting_places3[which(Mutasa5_clust_trips_shorter2_no_home_biting_places3$partid %in% Mutasa_parts_1_16),]
Mutasa5_clust_trips_shorter2_no_home_biting_places3_17_34 <- Mutasa5_clust_trips_shorter2_no_home_biting_places3[which(Mutasa5_clust_trips_shorter2_no_home_biting_places3$partid %in% Mutasa_parts_17_34),]
Mutasa5_clust_trips_shorter2_no_home_biting_places3_35_54 <- Mutasa5_clust_trips_shorter2_no_home_biting_places3[which(Mutasa5_clust_trips_shorter2_no_home_biting_places3$partid %in% Mutasa_parts_35_54),]
Mutasa5_clust_trips_shorter2_no_home_biting_places3_55_up <- Mutasa5_clust_trips_shorter2_no_home_biting_places3[which(Mutasa5_clust_trips_shorter2_no_home_biting_places3$partid %in% Mutasa_parts_55_up),]

favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
# 1  1      1  1  13 1.938 2.999 16       0
favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2  22 2.362 3.285 116       0
favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1 1.75  22 2.337 3.533 86       0
favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2  22 2.362 3.285 116       0

summary(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2))
#  1  2  3 13 
# 13  1  1  1 
summary(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2))
#   1  2  3  4  6  7  8  9 11 22 
#  72 20  8  4  3  1  3  2  1  2 
summary(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2))
#    1  2  3  4  5  7  8  9 10 11 14 22 
#   64  9  2  1  1  1  1  2  1  1  2  1 
summary(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2))
#  1  2  3  4 
# 26  6  2  1 

### PER PART ###
trip_count_median_1_16 <- favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_1_16$partid)$median
favstats(trip_count_median_1_16)
# min Q1 median Q3 max  mean  sd  n missing
#     1  1      1  1  13 2.444 3.972 9       0
trip_count_median_17_34 <- favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_17_34$partid)$median
favstats(trip_count_median_17_34)
# min Q1 median Q3 max  mean  sd  n missing
# 1  1    1.5  2   6 1.744 1.131 43       0
trip_count_median_35_54 <- favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_35_54$partid)$median
favstats(trip_count_median_35_54)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1 1.875  11    2 2.43 38       0
trip_count_median_55_up <- favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_55_up$partid)$median
favstats(trip_count_median_55_up)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  1   3 1.194 0.5185 18       0

trips_part2_age <- as.data.frame(cbind(c(trip_count_median_1_16, trip_count_median_17_34, trip_count_median_35_54, trip_count_median_55_up),
                                       c(rep("1-16",length(trip_count_median_1_16)),
                                         rep("17-34",length(trip_count_median_17_34)),
                                         rep("35-54",length(trip_count_median_35_54)),
                                         rep("55+",length(trip_count_median_55_up)))))
colnames(trips_part2_age) <- c("values","type")
trips_part2_age$values <- as.numeric(as.character(trips_part2_age$values))
trips_part2_age$type <- factor(trips_part2_age$type)

ggplot(data=trips_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,13.25)) + 
  scale_x_continuous(breaks=c(seq(0,14,2)))
#
##### percent biting time spent at home vs elsewhere #####
Mutasa5_clusts_shorter2_biting_places3b_1_16 <- Mutasa5_clusts_shorter2_biting_places3b[which(Mutasa5_clusts_shorter2_biting_places3b$partid %in% Mutasa_parts_1_16),]
Mutasa5_clusts_shorter2_biting_places3b_17_34 <- Mutasa5_clusts_shorter2_biting_places3b[which(Mutasa5_clusts_shorter2_biting_places3b$partid %in% Mutasa_parts_17_34),]
Mutasa5_clusts_shorter2_biting_places3b_35_54 <- Mutasa5_clusts_shorter2_biting_places3b[which(Mutasa5_clusts_shorter2_biting_places3b$partid %in% Mutasa_parts_35_54),]
Mutasa5_clusts_shorter2_biting_places3b_55_up <- Mutasa5_clusts_shorter2_biting_places3b[which(Mutasa5_clusts_shorter2_biting_places3b$partid %in% Mutasa_parts_55_up),]

favstats(Mutasa5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time)
favstats(Mutasa5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time[which(Mutasa5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time != 0)])
favstats(Mutasa5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time)
favstats(Mutasa5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time[which(Mutasa5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time != 0)])
favstats(Mutasa5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time)
favstats(Mutasa5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time[which(Mutasa5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time != 0)])
favstats(Mutasa5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time)
favstats(Mutasa5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time[which(Mutasa5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time != 0)])


perc_time_home2_age <- as.data.frame(cbind(c(Mutasa5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time, Mutasa5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time,
                                             Mutasa5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time, Mutasa5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time),
                                           c(rep("1-16",length(Mutasa5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time)),
                                             rep("17-34",length(Mutasa5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time)),
                                             rep("35-54",length(Mutasa5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time)),
                                             rep("55+",length(Mutasa5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time)))))
colnames(perc_time_home2_age) <- c("values","type")
perc_time_home2_age$values <- as.numeric(as.character(perc_time_home2_age$values))
perc_time_home2_age$type <- factor(perc_time_home2_age$type)

ggplot(data=perc_time_home2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-1,101))
##
###############
###############
##### NCHELENGE #####
###### General Age Analysis ####
Nchelenge_part_info <- Nchelenge5_clusts_shorter2[!duplicated(Nchelenge5_clusts_shorter2$partid),c(1:3)]
summary(Nchelenge_part_info$age_all)
# 10 - 72; mean=34.3, median=33, 1Q=20, 3Q=46.5
Nchelenge_part_info$age_18 <- factor(0, levels=c(0,1))
Nchelenge_part_info$age_18[which(Nchelenge_part_info$age_all >= 18)] <- 1
summary(as.factor(Nchelenge_part_info$age_18))
#  0   1
# 12  63  
Nchelenge_part_info$age_16_35_55_up <- factor(0, levels=c(0,1,2,3)) # 0 = <16; 1 = 17-34; 2 = 35-54; 3= 55+
Nchelenge_part_info$age_16_35_55_up[which(Nchelenge_part_info$age_all >= 17 & Nchelenge_part_info$age_all < 35)] <- 1
Nchelenge_part_info$age_16_35_55_up[which(Nchelenge_part_info$age_all >= 35 & Nchelenge_part_info$age_all < 55)] <- 2
Nchelenge_part_info$age_16_35_55_up[which(Nchelenge_part_info$age_all >= 55)] <- 3
summary(as.factor(Nchelenge_part_info$age_16_35_55_up))
#  0   1   2   3
# 11  28  27   9

ggplot(data=Nchelenge_part_info, aes(x=age_all)) +
  geom_histogram(binwidth = 1, aes(fill=age_16_35_55_up)) +
  coord_cartesian(xlim=c(0,75))
#
###### Split by Age ######
Nchelenge_parts_1_16 <- Nchelenge_part_info$partid[which(Nchelenge_part_info$age_16_35_55_up == 0)]
Nchelenge_parts_17_34 <- Nchelenge_part_info$partid[which(Nchelenge_part_info$age_16_35_55_up == 1)]
Nchelenge_parts_35_54 <- Nchelenge_part_info$partid[which(Nchelenge_part_info$age_16_35_55_up == 2)]
Nchelenge_parts_55_up <- Nchelenge_part_info$partid[which(Nchelenge_part_info$age_16_35_55_up == 3)]
##
### number of locations (without home) ####
Nchelenge5_clusts_shorter2_no_home_1_16 <- Nchelenge5_clusts_shorter2_no_home[which(Nchelenge5_clusts_shorter2_no_home$partid %in% Nchelenge_parts_1_16),]
Nchelenge5_clusts_shorter2_no_home_17_34 <- Nchelenge5_clusts_shorter2_no_home[which(Nchelenge5_clusts_shorter2_no_home$partid %in% Nchelenge_parts_17_34),]
Nchelenge5_clusts_shorter2_no_home_35_54 <- Nchelenge5_clusts_shorter2_no_home[which(Nchelenge5_clusts_shorter2_no_home$partid %in% Nchelenge_parts_35_54),]
Nchelenge5_clusts_shorter2_no_home_55_up <- Nchelenge5_clusts_shorter2_no_home[which(Nchelenge5_clusts_shorter2_no_home$partid %in% Nchelenge_parts_55_up),]

loc_counts_1_16 <- favstats(Nchelenge5_clusts_shorter2_no_home_1_16$clust~Nchelenge5_clusts_shorter2_no_home_1_16$partid)$n
only_home_1_16 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_1_16)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_1_16$partid))
loc_counts_1_162 <- c(loc_counts_1_16,rep(0, only_home_1_16))
favstats(loc_counts_1_162)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0  2      5 8.5  11 5.090909 3.753786 11       0
loc_counts_17_34 <- favstats(Nchelenge5_clusts_shorter2_no_home_17_34$clust~Nchelenge5_clusts_shorter2_no_home_17_34$partid)$n
only_home_17_34 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_17_34)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_17_34$partid))
loc_counts_17_342 <- c(loc_counts_17_34,rep(0, only_home_17_34))
favstats(loc_counts_17_342)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0 4.75    8.5 14.25  23 9.357143 5.970385 28       0
loc_counts_35_54 <- favstats(Nchelenge5_clusts_shorter2_no_home_35_54$clust~Nchelenge5_clusts_shorter2_no_home_35_54$partid)$n
only_home_35_54 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_35_54)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_35_54$partid))
loc_counts_35_542 <- c(loc_counts_35_54,rep(0, only_home_35_54))
favstats(loc_counts_35_542)
# min   Q1 median     Q3  max     mean        sd  n missing
#   0 5.5      8 12  19 9.074074 5.413132 27       0
loc_counts_55_up <- favstats(Nchelenge5_clusts_shorter2_no_home_55_up$clust~Nchelenge5_clusts_shorter2_no_home_55_up$partid)$n
only_home_55_up <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_55_up)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_55_up$partid))
loc_counts_55_up2 <- c(loc_counts_55_up,rep(0, only_home_55_up))
favstats(loc_counts_55_up2)
# min   Q1 median     Q3  max     mean        sd  n missing
#   0  3     11 14  19 8.888889 6.808899 9       0

locs_count2 <- as.data.frame(cbind(c(loc_counts_1_162, loc_counts_17_342, loc_counts_35_542, loc_counts_55_up2),
                                   c(rep("1-16",length(loc_counts_1_162)),
                                     rep("17-34",length(loc_counts_17_342)),
                                     rep("35-54",length(loc_counts_35_542)),
                                     rep("55+",length(loc_counts_55_up2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

shapiro.test(locs_count2$values) ## not normal
kruskal.test(locs_count2$values, locs_count2$type) # p = 0.18
pairwise.wilcox.test(locs_count2$values, locs_count2$type, p.adjust.method = "holm")

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-2,26)) +
  scale_x_continuous(breaks=c(seq(0,25,5)))
#
### distance of locations (without home) ####
### TOTAL ###
favstats(Nchelenge5_clusts_shorter2_no_home_1_16$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
# 0.03850304 0.1946278 0.5277768 5.273784 42.07533 5.511531 10.12079 56       0
favstats(Nchelenge5_clusts_shorter2_no_home_17_34$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
#  0 0.3698045 1.175084 4.180639 146.5504 8.944337 25.45071 262       0
favstats(Nchelenge5_clusts_shorter2_no_home_35_54$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
# 0.009427688 0.395577 1.487539 5.454201 187.2311 6.166172 18.28843 245       0
favstats(Nchelenge5_clusts_shorter2_no_home_55_up$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 0.02460256 0.3555845 1.084418 3.20567 20.27392 3.010719 4.328903 80       0

### PER PART ###
medians_km_1_16 <- (favstats(Nchelenge5_clusts_shorter2_no_home_1_16$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_1_16$partid)$median)
medians_km_17_34 <- (favstats(Nchelenge5_clusts_shorter2_no_home_17_34$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_17_34$partid)$median)
medians_km_35_54 <- (favstats(Nchelenge5_clusts_shorter2_no_home_35_54$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_35_54$partid)$median)
medians_km_55_up <- (favstats(Nchelenge5_clusts_shorter2_no_home_55_up$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_55_up$partid)$median)

favstats(medians_km_1_16) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.1811803 0.2862432 0.5790771 1.646362 25.71711 3.633996 7.913153 10       0
favstats(medians_km_17_34) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
# 0.1012507 0.6154615 1.166392 2.614608 85.42321 5.006667 16.23232 27       0
favstats(medians_km_35_54) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.1602008 0.8441343 1.489453 4.293316 186.7748 11.14915 36.32534 26       0
favstats(medians_km_55_up) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
# 0.3700853 0.5344191 0.7619818 1.571078 3.167973 1.200207 0.9905214 8       0

## REMOVE OUTLIER ##
medians_km_17_34_short <- medians_km_17_34[-which(medians_km_17_34 > 70)]
favstats(medians_km_17_34_short)
# min        Q1    median       Q3      max     mean       sd  n missing
# 0.1012507 0.6149524 1.160551 2.178421 11.50774 1.913723 2.325044 26       0

medians_km_35_54_short <- medians_km_35_54[-which(medians_km_35_54 > 70)]
favstats(medians_km_35_54_short)
# min        Q1    median       Q3      max     mean       sd  n missing
# 0.1602008 0.8181657 1.312188 3.403919 24.49658 4.124123 6.158095 25       0

dists_km2_age <- as.data.frame(cbind(c(medians_km_1_16, medians_km_17_34_short, medians_km_35_54_short, medians_km_55_up),
                                     c(rep("1-16",length(medians_km_1_16)),
                                       rep("17-34",length(medians_km_17_34_short)),
                                       rep("35-54",length(medians_km_35_54_short)),
                                       rep("55+",length(medians_km_55_up)))))
colnames(dists_km2_age) <- c("values","type")
dists_km2_age$values <- as.numeric(as.character(dists_km2_age$values))
dists_km2_age$type <- factor(dists_km2_age$type)

shapiro.test(dists_km2_age$values) ## not normal
kruskal.test(dists_km2_age$values, dists_km2_age$type) # p = 0.35
pairwise.wilcox.test(dists_km2_age$values, dists_km2_age$type, p.adjust.method = "holm")

ggplot(data=dists_km2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-1,30))
#
###
### number of trips (without home) ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2_no_home3_1_16 <- Nchelenge5_clust_trips_shorter2_no_home3[which(Nchelenge5_clust_trips_shorter2_no_home3$partid %in% Nchelenge_parts_1_16),]
Nchelenge5_clust_trips_shorter2_no_home3_17_34 <- Nchelenge5_clust_trips_shorter2_no_home3[which(Nchelenge5_clust_trips_shorter2_no_home3$partid %in% Nchelenge_parts_17_34),]
Nchelenge5_clust_trips_shorter2_no_home3_35_54 <- Nchelenge5_clust_trips_shorter2_no_home3[which(Nchelenge5_clust_trips_shorter2_no_home3$partid %in% Nchelenge_parts_35_54),]
Nchelenge5_clust_trips_shorter2_no_home3_55_up <- Nchelenge5_clust_trips_shorter2_no_home3[which(Nchelenge5_clust_trips_shorter2_no_home3$partid %in% Nchelenge_parts_55_up),]

favstats(Nchelenge5_clust_trips_shorter2_no_home3_1_16$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
# 1  1      2 2.25   9 2.089286 1.621187 56       0
favstats(Nchelenge5_clust_trips_shorter2_no_home3_17_34$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      2  3  23 2.798479 2.966942 263       0
favstats(Nchelenge5_clust_trips_shorter2_no_home3_35_54$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      1  2  20 2.432653 2.898656 245       0
favstats(Nchelenge5_clust_trips_shorter2_no_home3_55_up$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      2  3   8 2.1375 1.5323 80       0

### PER PART ###
trips_part_median_1_16 <- favstats(Nchelenge5_clust_trips_shorter2_no_home3_1_16$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_1_16$partid)$median
trips_part_median_17_34 <- favstats(Nchelenge5_clust_trips_shorter2_no_home3_17_34$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_17_34$partid)$median
trips_part_median_35_54 <- favstats(Nchelenge5_clust_trips_shorter2_no_home3_35_54$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_35_54$partid)$median
trips_part_median_55_up <- favstats(Nchelenge5_clust_trips_shorter2_no_home3_55_up$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_55_up$partid)$median

favstats(trips_part_median_1_16)
#   min  Q1   median   Q3      max       mean       sd  n missing
#    1 1.125      2  2   3 1.75 0.6346478 10       0
favstats(trips_part_median_17_34)
#   min Q1 median  Q3   max      mean       sd  n missing
#  1  1      2  2   4 1.740741 0.6703956 27       0
favstats(trips_part_median_35_54)
#   min  Q1   median   Q3      max       mean       sd  n missing
#  1  1    1.5  2   3 1.538462 0.5817745 26       0
favstats(trips_part_median_55_up)
#   min Q1 median  Q3   max      mean       sd  n missing
#  1  1    1.5  2   2  1.5 0.5345225 8       0

trips_part2_age <- as.data.frame(cbind(c(trips_part_median_1_16, trips_part_median_17_34, trips_part_median_35_54, trips_part_median_55_up),
                                       c(rep("1-16",length(trips_part_median_1_16)),
                                         rep("17-34",length(trips_part_median_17_34)),
                                         rep("35-54",length(trips_part_median_35_54)),
                                         rep("55+",length(trips_part_median_55_up)))))
colnames(trips_part2_age) <- c("values","type")
trips_part2_age$values <- as.numeric(as.character(trips_part2_age$values))
trips_part2_age$type <- factor(trips_part2_age$type)
kruskal.test(trips_part2_age$values, trips_part2_age$type) # p = 0.35
pairwise.wilcox.test(trips_part2_age$values, trips_part2_age$type, p.adjust.method = "holm")

ggplot(data=trips_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,6)) +
  scale_x_continuous(breaks=c(0:6))
#
### time per trip ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2_no_home_1_16 <- Nchelenge5_clust_trips_shorter2_no_home[which(Nchelenge5_clust_trips_shorter2_no_home$partid %in% Nchelenge_parts_1_16),]
Nchelenge5_clust_trips_shorter2_no_home_17_34 <- Nchelenge5_clust_trips_shorter2_no_home[which(Nchelenge5_clust_trips_shorter2_no_home$partid %in% Nchelenge_parts_17_34),]
Nchelenge5_clust_trips_shorter2_no_home_35_54 <- Nchelenge5_clust_trips_shorter2_no_home[which(Nchelenge5_clust_trips_shorter2_no_home$partid %in% Nchelenge_parts_35_54),]
Nchelenge5_clust_trips_shorter2_no_home_55_up <- Nchelenge5_clust_trips_shorter2_no_home[which(Nchelenge5_clust_trips_shorter2_no_home$partid %in% Nchelenge_parts_55_up),]

favstats(Nchelenge5_clust_trips_shorter2_no_home_1_16$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2558333 0.6358333 1.444444 3.673889 329.3114 10.08758 40.08207 117       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_17_34$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2416667 0.5775 1.229167 2.753333 379.3781 5.867487 22.89963 727       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_35_54$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2425 0.57375 1.143333 2.575694 394.7653 5.418526 24.09383 595       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_55_up$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2441667 0.5456944 1.143889 2.374306 122.0117 2.773502 9.723764 171       0

### PER PART ###
trip_time_part_median_1_16 <- favstats(Nchelenge5_clust_trips_shorter2_no_home_1_16$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_1_16$partid)$median
trip_time_part_median_17_34 <- favstats(Nchelenge5_clust_trips_shorter2_no_home_17_34$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_17_34$partid)$median
trip_time_part_median_35_54 <- favstats(Nchelenge5_clust_trips_shorter2_no_home_35_54$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_35_54$partid)$median
trip_time_part_median_55_up <- favstats(Nchelenge5_clust_trips_shorter2_no_home_55_up$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_55_up$partid)$median

favstats(trip_time_part_median_1_16)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4844444 0.9115972 1.768889 2.157083 2.747361 1.642486 0.7951069 10       0
favstats(trip_time_part_median_17_34)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.4966667 0.8788889 1.382917 1.720625 6.663194 1.581893 1.189569 27       0
favstats(trip_time_part_median_35_54)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4955556 0.8959375 1.170694 1.71184 6.148333 1.562196 1.21398 26       0
favstats(trip_time_part_median_55_up)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.7569444 1.012014 1.25125 1.607882 3.828056 1.547587 0.9833209 8       0

trip_times_part2_age <- as.data.frame(cbind(c(trip_time_part_median_1_16, trip_time_part_median_17_34, trip_time_part_median_35_54, trip_time_part_median_55_up),
                                            c(rep("1-16",length(trip_time_part_median_1_16)),
                                              rep("17-34",length(trip_time_part_median_17_34)),
                                              rep("35-54",length(trip_time_part_median_35_54)),
                                              rep("55+",length(trip_time_part_median_55_up)))))
colnames(trip_times_part2_age) <- c("values","type")
trip_times_part2_age$values <- as.numeric(as.character(trip_times_part2_age$values))
trip_times_part2_age$type <- factor(trip_times_part2_age$type)

kruskal.test(trip_times_part2_age$values, trip_times_part2_age$type) # p = 0.35
pairwise.wilcox.test(trip_times_part2_age$values, trip_times_part2_age$type, p.adjust.method = "holm")

ggplot(data=trip_times_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge",  binwidth = 0.5, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,7))+
  scale_x_continuous(breaks=c(0:7))
#
### Home Loc ####
Nchelenge5_clusts_shorter2_home_1_16 <- Nchelenge5_clusts_shorter2_home[which(Nchelenge5_clusts_shorter2_home$partid %in% Nchelenge_parts_1_16),] 
no_home_1_16 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_1_16)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_home_1_16$partid))# 1
## none with no home
Nchelenge5_clusts_shorter2_home_17_34 <- Nchelenge5_clusts_shorter2_home[which(Nchelenge5_clusts_shorter2_home$partid %in% Nchelenge_parts_17_34),] 
no_home_17_34 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_17_34)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_home_17_34$partid))# 1
## none with no home
Nchelenge5_clusts_shorter2_home_35_54 <- Nchelenge5_clusts_shorter2_home[which(Nchelenge5_clusts_shorter2_home$partid %in% Nchelenge_parts_35_54),] 
no_home_35_54 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_35_54)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_home_35_54$partid))# 1
## one with no home
Nchelenge5_clusts_shorter2_home_55_up <- Nchelenge5_clusts_shorter2_home[which(Nchelenge5_clusts_shorter2_home$partid %in% Nchelenge_parts_55_up),] 
no_home_55_up <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_55_up)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_home_55_up$partid))# 1
## none with no home

### number of 'trips' to home loc ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2_home3_1_16 <- Nchelenge5_clust_trips_shorter2_home3[which(Nchelenge5_clust_trips_shorter2_home3$partid %in% Nchelenge_parts_1_16),]
Nchelenge5_clust_trips_shorter2_home3_17_34 <- Nchelenge5_clust_trips_shorter2_home3[which(Nchelenge5_clust_trips_shorter2_home3$partid %in% Nchelenge_parts_17_34),]
Nchelenge5_clust_trips_shorter2_home3_35_54 <- Nchelenge5_clust_trips_shorter2_home3[which(Nchelenge5_clust_trips_shorter2_home3$partid %in% Nchelenge_parts_35_54),]
Nchelenge5_clust_trips_shorter2_home3_55_up <- Nchelenge5_clust_trips_shorter2_home3[which(Nchelenge5_clust_trips_shorter2_home3$partid %in% Nchelenge_parts_55_up),]

favstats(Nchelenge5_clust_trips_shorter2_home3_1_16$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#    1  2      3 8.5  15 5.636364 4.945154 11       0
favstats(Nchelenge5_clust_trips_shorter2_home3_17_34$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#      1 3.75      7 14.5  48 11.10714 11.24258 28       0
favstats(Nchelenge5_clust_trips_shorter2_home3_35_54$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#    1 4.5      9 13.5  34 10.85185 8.834825 27       0
favstats(Nchelenge5_clust_trips_shorter2_home3_55_up$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#    1  5     17 19  30 13.55556 9.72254 9       0
#
### total time at home loc ####
Nchelenge5_clust_trips_shorter2_home4_1_16 <- Nchelenge5_clust_trips_shorter2_home4[which(Nchelenge5_clust_trips_shorter2_home4$partid %in% Nchelenge_parts_1_16),]
favstats(Nchelenge5_clust_trips_shorter2_home4_1_16$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
#  0.260636 12.74333 22.18398 26.38974 31.1795 19.21422 10.34955 11       0
Nchelenge5_clust_trips_shorter2_home4_17_34 <- Nchelenge5_clust_trips_shorter2_home4[which(Nchelenge5_clust_trips_shorter2_home4$partid %in% Nchelenge_parts_17_34),]
favstats(Nchelenge5_clust_trips_shorter2_home4_17_34$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
#   0.028894 12.54172 22.86267 28.61703 40.73759 20.64817 10.97598 28       0
Nchelenge5_clust_trips_shorter2_home4_35_54 <- Nchelenge5_clust_trips_shorter2_home4[which(Nchelenge5_clust_trips_shorter2_home4$partid %in% Nchelenge_parts_35_54),]
favstats(Nchelenge5_clust_trips_shorter2_home4_35_54$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
# 1.588992 19.58215 24.6938 31.41483 36.04353 23.49449 9.446279 27       0
Nchelenge5_clust_trips_shorter2_home4_55_up <- Nchelenge5_clust_trips_shorter2_home4[which(Nchelenge5_clust_trips_shorter2_home4$partid %in% Nchelenge_parts_55_up),]
favstats(Nchelenge5_clust_trips_shorter2_home4_55_up$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
# 9.579859 16.46552 24.06558 24.49137 37.21134 22.86481 8.489825 9       0
#
### percent time at home ####
favstats(Nchelenge5_clusts_shorter2_home_1_16$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
#  11.62073 87.4981 97.89361 99.36384 100 87.10382 25.90311 11       0
favstats(Nchelenge5_clusts_shorter2_home_17_34$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#  1.618541 70.31373 91.75945 94.98805 100 80.54125 22.95651 27       0
favstats(Nchelenge5_clusts_shorter2_home_35_54$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
#  12.92136 80.469 90.16104 97.85489 100 83.99721 19.86029 27       0
favstats(Nchelenge5_clusts_shorter2_home_55_up$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#  81.23371 87.1791 91.73057 98.12107 100 91.88917 7.090506 9       0

perc_home2_age <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_home_1_16$perc_time_home, Nchelenge5_clusts_shorter2_home_17_34$perc_time_home, Nchelenge5_clusts_shorter2_home_35_54$perc_time_home, Nchelenge5_clusts_shorter2_home_55_up$perc_time_home),
                                      c(rep("1-16",length(Nchelenge5_clusts_shorter2_home_1_16$perc_time_home)),
                                        rep("17-34",length(Nchelenge5_clusts_shorter2_home_17_34$perc_time_home)),
                                        rep("35-54",length(Nchelenge5_clusts_shorter2_home_35_54$perc_time_home)),
                                        rep("55+",length(Nchelenge5_clusts_shorter2_home_55_up$perc_time_home)))))
colnames(perc_home2_age) <- c("values","type")
perc_home2_age$values <- as.numeric(as.character(perc_home2_age$values))
perc_home2_age$type <- factor(perc_home2_age$type)
kruskal.test(perc_home2_age$values, perc_home2_age$type) # p = 0.35
pairwise.wilcox.test(perc_home2_age$values, perc_home2_age$type, p.adjust.method = "holm")

ggplot(data=perc_home2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,101))
##
### Average Direction of Travel #####
temp <- Nchelenge5_clusts_shorter2_no_home[!duplicated(Nchelenge5_clusts_shorter2_no_home$partid),]

age_group <- Nchelenge5_clusts_shorter2_no_home$age_all[!duplicated(Nchelenge5_clusts_shorter2_no_home$partid)]

temp$age_group <- factor(0, levels=c(0:3)) # 0 = <16; 1 = 17-34; 2 = 35-54; 3= 55+
temp$age_group[which(age_group >= 17 & age_group < 35)] <- 1
temp$age_group[which(age_group >= 35 & age_group < 55)] <- 2
temp$age_group[which(age_group >= 55)] <- 3
summary(as.factor(temp$age_group))

ggplot(data=temp, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=age_group), position=position_dodge2(preserve = "single"), width=6) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Age Group", values=agecols[c(5,4,2,1)], labels=c("55+", "35-54", "17-34", "0-16"))
#
favstats(temp$avg_dir_from_home_bearing~temp$age_group)
# temp$age_group       min        Q1   median       Q3      max     mean       sd  n missing
# 1              3 48.434751  89.08409 157.8047 224.3082 285.5163 161.5314 84.85274  8       0
# 2              2 15.655534 127.93779 200.7787 284.4924 358.9609 203.0965 91.08068 26       0
# 3              1  7.761391 137.89064 208.6944 305.6945 347.7768 216.8442 99.15089 27       0
# 4              0 88.800230 147.19618 227.5121 273.2939 357.0854 217.3343 93.29801 10       0
favstats(temp$avg_dir_from_home_mag~temp$age_group)
# temp$age_group       min        Q1    median       Q3        max     mean        sd  n
# 1              3 0.4501813 0.5330575 0.7348198 2.642140   3.580177 1.524907  1.334963  8
# 2              2 0.1602008 0.7095076 1.8278507 7.023395 105.150857 8.456090 20.579861 26
# 3              1 0.1034942 0.6987528 1.4386150 2.667428  74.437652 5.330296 14.297053 27
# 4              0 0.1613577 0.2671998 1.1108374 4.319618  21.569291 4.932807  7.654216 10

temp2 <- temp[which(temp$avg_dir_from_home_mag < 65),] ## 164 of 176 parts

favstats(temp2$avg_dir_from_home_bearing~temp2$age_group)
#   temp2$age_group       min        Q1   median       Q3      max     mean       sd  n
# 1               3 48.434751  89.08409 157.8047 224.3082 285.5163 161.5314 84.85274  8
# 2               2 15.655534 120.04147 217.2187 289.4638 358.9609 204.0885 92.81537 25
# 3               1  7.761391 137.63183 221.6625 306.5534 347.7768 219.8437 99.85740 26
# 4               0 88.800230 147.19618 227.5121 273.2939 357.0854 217.3343 93.29801 10
favstats(temp2$avg_dir_from_home_mag~temp2$age_group)
#   temp2$age_group       min        Q1    median       Q3       max     mean       sd  n
# 1               3 0.4501813 0.5330575 0.7348198 2.642140  3.580177 1.524907 1.334963  8
# 2               2 0.1602008 0.6688438 1.7698629 6.130227 20.624539 4.588299 6.001366 25
# 3               1 0.1034942 0.6752155 1.3598931 2.410679 16.646568 2.672320 3.768480 26
# 4               0 0.1613577 0.2671998 1.1108374 4.319618 21.569291 4.932807 7.654216 10

kruskal.test(temp2$avg_dir_from_home_bearing, temp2$age_group) # p = 0.35
pairwise.wilcox.test(temp2$avg_dir_from_home_bearing, temp2$age_group, p.adjust.method = "holm")

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=age_group), position=position_dodge2(preserve = "single"), width=0.5) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Age Group", values=agecols[c(5,4,2,1)], labels=c("55+", "35-54", "17-34", "0-16"))




temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)] <- -(360-temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)])
ggplot(data=temp, aes(x=avg_dir_from_home_values, y=..count.., fill=age_group)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~age_group, labeller = labeller(age_group = c("0"="0-16", "1"="17-34","2"="35-54", "3"="55+")))+
  scale_fill_manual(values=agecols[c(1,2,4,5)]) +
  scale_y_continuous(breaks=c(seq(0,8,2)))


#
### time spent at clusters vs elsewhere (travel, less important places, etc) ####
# percent time in clusters
Nchelenge5_clusts_shorter2_1_16 <- Nchelenge5_clusts_shorter2[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_1_16),]
Nchelenge5_clusts_shorter2_17_34 <- Nchelenge5_clusts_shorter2[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_17_34),]
Nchelenge5_clusts_shorter2_35_54 <- Nchelenge5_clusts_shorter2[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_35_54),]
Nchelenge5_clusts_shorter2_55_up <- Nchelenge5_clusts_shorter2[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_55_up),]

Nchelenge5_clusts_shorter2_1_16$perc_time_in_all_clusts <- Nchelenge5_clusts_shorter2_1_16$part_trips_time_new/Nchelenge5_clusts_shorter2_1_16$total_time
Nchelenge5_clusts_shorter2_17_34$perc_time_in_all_clusts <- Nchelenge5_clusts_shorter2_17_34$part_trips_time_new/Nchelenge5_clusts_shorter2_17_34$total_time
Nchelenge5_clusts_shorter2_35_54$perc_time_in_all_clusts <- Nchelenge5_clusts_shorter2_35_54$part_trips_time_new/Nchelenge5_clusts_shorter2_35_54$total_time
Nchelenge5_clusts_shorter2_55_up$perc_time_in_all_clusts <- Nchelenge5_clusts_shorter2_55_up$part_trips_time_new/Nchelenge5_clusts_shorter2_55_up$total_time

summary(Nchelenge5_clusts_shorter2_1_16$perc_time_in_all_clusts[!duplicated(Nchelenge5_clusts_shorter2_1_16$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8451  0.8610  0.9074  0.9118  0.9541  0.9932 
summary(Nchelenge5_clusts_shorter2_17_34$perc_time_in_all_clusts[!duplicated(Nchelenge5_clusts_shorter2_17_34$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8663  0.9148  0.9509  0.9451  0.9895  0.9955 
summary(Nchelenge5_clusts_shorter2_35_54$perc_time_in_all_clusts[!duplicated(Nchelenge5_clusts_shorter2_35_54$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6629  0.8749  0.9248  0.9012  0.9642  0.9984 
summary(Nchelenge5_clusts_shorter2_55_up$perc_time_in_all_clusts[!duplicated(Nchelenge5_clusts_shorter2_55_up$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7781  0.8729  0.9317  0.9195  0.9782  0.9950 


# time outside clusters
summary(Nchelenge5_clusts_shorter2_1_16$total_time[!duplicated(Nchelenge5_clusts_shorter2_1_16$partid)] - Nchelenge5_clusts_shorter2_1_16$part_trips_time_new[!duplicated(Nchelenge5_clusts_shorter2_1_16$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2509  1.2039  2.9537  2.6542  4.3409  4.8027 
summary(Nchelenge5_clusts_shorter2_17_34$total_time[!duplicated(Nchelenge5_clusts_shorter2_17_34$partid)] - Nchelenge5_clusts_shorter2_17_34$part_trips_time_new[!duplicated(Nchelenge5_clusts_shorter2_17_34$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1395  0.4087  1.4091  1.8793  3.0862  4.7086 
summary(Nchelenge5_clusts_shorter2_35_54$total_time[!duplicated(Nchelenge5_clusts_shorter2_35_54$partid)] - Nchelenge5_clusts_shorter2_35_54$part_trips_time_new[!duplicated(Nchelenge5_clusts_shorter2_35_54$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2275  1.1145  2.0146  3.3552  4.5656  9.8070 
summary(Nchelenge5_clusts_shorter2_55_up$total_time[!duplicated(Nchelenge5_clusts_shorter2_55_up$partid)] - Nchelenge5_clusts_shorter2_55_up$part_trips_time_new[!duplicated(Nchelenge5_clusts_shorter2_55_up$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2383  0.7049  2.3483  2.7655  4.5573  7.0787 
#
### Time at locations ####
favstats(Nchelenge5_clusts_shorter2_no_home_1_16$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.4011111 1.320764 2.494028 6.248264 630.46 21.07584 85.23548 56       0
favstats(Nchelenge5_clusts_shorter2_no_home_17_34$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2480556 1.157153 2.518611 7.09375 494.26 16.26116 56.64596 262       0
favstats(Nchelenge5_clusts_shorter2_no_home_35_54$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2425 0.9313889 1.797778 4.078611 784.0228 13.11252 60.65909 245       0
favstats(Nchelenge5_clusts_shorter2_no_home_55_up$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2583333 0.8458333 2.325694 5.476875 124.6425 5.928361 14.90274 80       0

### PER PART ###
clust_time_part_median_1_16 <- (favstats(Nchelenge5_clusts_shorter2_no_home_1_16$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_1_16$partid)$median)*24
clust_time_part_median_17_34 <- (favstats(Nchelenge5_clusts_shorter2_no_home_17_34$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_17_34$partid)$median)*24
clust_time_part_median_35_54 <- (favstats(Nchelenge5_clusts_shorter2_no_home_35_54$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_35_54$partid)$median)*24
clust_time_part_median_55_up <- (favstats(Nchelenge5_clusts_shorter2_no_home_55_up$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_55_up$partid)$median)*24

favstats(clust_time_part_median_1_16)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 1.345833 1.794757 2.781111 3.231319 4.030278 2.612236 0.9289645 10       0
favstats(clust_time_part_median_17_34)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.8686111 1.685278 2.4625 3.534861 12.90486 3.187757 2.525574 27       0
favstats(clust_time_part_median_35_54)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.9338889 1.496389 1.736181 3.000174 17.23194 3.033168 3.620338 26       0
favstats(clust_time_part_median_55_up)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.4605556 2.0225 2.572292 3.10625 3.828056 2.472135 1.030378 8       0


clust_time_part_age <- as.data.frame(cbind(c(clust_time_part_median_1_16, clust_time_part_median_17_34,clust_time_part_median_35_54,clust_time_part_median_55_up),
                                           c(rep("1-16",length(clust_time_part_median_1_16)),
                                             rep("17-34",length(clust_time_part_median_17_34)),
                                             rep("35-54",length(clust_time_part_median_35_54)),
                                             rep("55+",length(clust_time_part_median_55_up)))))
colnames(clust_time_part_age) <- c("values","type")
clust_time_part_age$values <- as.numeric(as.character(clust_time_part_age$values))
clust_time_part_age$type <- factor(clust_time_part_age$type)
kruskal.test(clust_time_part_age$values, clust_time_part_age$type) # p = 0.35
pairwise.wilcox.test(clust_time_part_age$values, clust_time_part_age$type, p.adjust.method = "holm")

ggplot(data=clust_time_part_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,20)) +
  scale_x_continuous(breaks=c(seq(0,20,2)))
#
#
#
####### Biting time metrics #####
### number locations outside home ####
Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16 <- Nchelenge5_clusts_shorter2_no_home_biting_places2[which(Nchelenge5_clusts_shorter2_no_home_biting_places2$partid %in% Nchelenge_parts_1_16),]
Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34 <- Nchelenge5_clusts_shorter2_no_home_biting_places2[which(Nchelenge5_clusts_shorter2_no_home_biting_places2$partid %in% Nchelenge_parts_17_34),]
Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54 <- Nchelenge5_clusts_shorter2_no_home_biting_places2[which(Nchelenge5_clusts_shorter2_no_home_biting_places2$partid %in% Nchelenge_parts_35_54),]
Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up <- Nchelenge5_clusts_shorter2_no_home_biting_places2[which(Nchelenge5_clusts_shorter2_no_home_biting_places2$partid %in% Nchelenge_parts_55_up),]

loc_counts_1_16 <- favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16$partid)$n
loc_counts_17_34 <- favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34$partid)$n
loc_counts_35_54 <- favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54$partid)$n
loc_counts_55_up <- favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up$partid)$n
only_home_1_16 <- nlevels(as.factor(Nchelenge5_clusts_shorter2_1_16$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16$partid))
only_home_17_34 <- nlevels(as.factor(Nchelenge5_clusts_shorter2_17_34$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34$partid))
only_home_35_54 <- nlevels(as.factor(Nchelenge5_clusts_shorter2_35_54$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54$partid))
only_home_55_up <- nlevels(as.factor(Nchelenge5_clusts_shorter2_55_up$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up$partid))
loc_counts_1_162 <- c(loc_counts_1_16,rep(0, only_home_1_16))
loc_counts_17_342 <- c(loc_counts_17_34,rep(0, only_home_17_34))
loc_counts_35_542 <- c(loc_counts_35_54,rep(0, only_home_35_54))
loc_counts_55_up2 <- c(loc_counts_55_up,rep(0, only_home_55_up))

favstats(loc_counts_1_162)
# min Q1 median Q3 max     mean       sd  n missing
#    0  0      0 2.5   6 1.636 2.378 11       0
favstats(loc_counts_17_342)
# min Q1 median Q3 max     mean       sd  n missing
#   0 0.75      1 3.25  12 2.429 2.974 28       0
favstats(loc_counts_35_542)
# min Q1 median Q3 max     mean       sd  n missing
#    0  0      1  2   6 1.333 1.468 27       0
favstats(loc_counts_55_up2)
# min Q1 median Q3 max     mean       sd  n missing
#    0  0      0  1   4 0.7778 1.302 9       0

summary(as.factor(loc_counts_1_162))
# 0 1 2 3 6 
# 6 1 1 1 2
summary(as.factor(loc_counts_17_342))
#  0  1  2  3  4  5  8  9 12 
# 7  8  4  2  3  1  1  1  1 
summary(as.factor(loc_counts_35_542))
# 0  1  2  3  4  6 
# 10  6  7  2  1  1 
summary(as.factor(loc_counts_55_up2))
# 0 1 4 
# 5 3 1 
## all those with 0 have biting time at home location (not just nowhere)

locs_count2 <- as.data.frame(cbind(c(loc_counts_1_162, loc_counts_17_342, loc_counts_35_542, loc_counts_55_up2),
                                   c(rep("1-16",length(loc_counts_1_162)),
                                     rep("17-34",length(loc_counts_17_342)),
                                     rep("35-54",length(loc_counts_35_542)),
                                     rep("55+",length(loc_counts_55_up2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) + 
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-0.65,12.35)) +
  scale_x_continuous(breaks=c(seq(0,12,1)))
#
##
### time per location ####
favstats((Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.9167 1.412  8.761 19.74 181.7 21.14 41.87 18       0
favstats((Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5112 1.731  5.711 16.69 193.3 19.32 36.04 68       0
favstats((Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5714 1.563  7.754 32.63 285.4 24.9 50.13 36       0
favstats((Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.7069 1.297  9.063 13.47 33.89 10.46 12.02 7       0

clusts_median_Nchelenge_1_16 <- (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16$partid)$median)*24
clusts_median_Nchelenge_17_34 <- (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34$partid)$median)*24
clusts_median_Nchelenge_35_54 <- (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54$partid)$median)*24
clusts_median_Nchelenge_55_up <- (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up$partid)$median)*24
favstats(clusts_median_Nchelenge_1_16)
# min       Q1   median  Q3      max     mean       sd  n missing
# 1.189 1.71   7.97 12.23 26.99 10.02 10.54 5       0
favstats(clusts_median_Nchelenge_17_34)
# min       Q1   median  Q3      max     mean       sd  n missing
# 1.066 3.189  8.167 17.27 82.97 14.42 18.88 21       0
favstats(clusts_median_Nchelenge_35_54)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.8406 7.973  11.42 34.98 79.21 23.08 23.05 17       0
favstats(clusts_median_Nchelenge_55_up)
# min       Q1   median  Q3      max     mean       sd  n missing
#  0.8618 1.514  5.452 15.35 33.89 11.41 15.44 4       0

clust_time_part_age <- as.data.frame(cbind(c(clusts_median_Nchelenge_1_16, clusts_median_Nchelenge_17_34, clusts_median_Nchelenge_35_54, clusts_median_Nchelenge_55_up),
                                           c(rep("1-16",length(clusts_median_Nchelenge_1_16)),
                                             rep("17-34",length(clusts_median_Nchelenge_17_34)),
                                             rep("35-54",length(clusts_median_Nchelenge_35_54)),
                                             rep("55+",length(clusts_median_Nchelenge_55_up)))))
colnames(clust_time_part_age) <- c("values","type")
clust_time_part_age$values <- as.numeric(as.character(clust_time_part_age$values))
clust_time_part_age$type <- factor(clust_time_part_age$type)

ggplot(data=clust_time_part_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-10,100)) +
  scale_x_continuous(breaks=c(seq(0,100,20)))
#
### number trips ####
Nchelenge5_clust_trips_shorter2_no_home_biting_places3_1_16 <- Nchelenge5_clust_trips_shorter2_no_home_biting_places3[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$partid %in% Nchelenge_parts_1_16),]
Nchelenge5_clust_trips_shorter2_no_home_biting_places3_17_34 <- Nchelenge5_clust_trips_shorter2_no_home_biting_places3[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$partid %in% Nchelenge_parts_17_34),]
Nchelenge5_clust_trips_shorter2_no_home_biting_places3_35_54 <- Nchelenge5_clust_trips_shorter2_no_home_biting_places3[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$partid %in% Nchelenge_parts_35_54),]
Nchelenge5_clust_trips_shorter2_no_home_biting_places3_55_up <- Nchelenge5_clust_trips_shorter2_no_home_biting_places3[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$partid %in% Nchelenge_parts_55_up),]

favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#    1  1      1 1.75   6 1.667 1.372 18       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#    1  1      1  2  13 2.182 2.517 66       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      2  3   9 2.371 2.102 35       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  1   2 1.143 0.378 7       0

summary(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2))
#  1  2  3  4  6 
# 13  2  1  1  1 
summary(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2))
#   1  2  3  4  5  6  9 10 12 13 
# 41 11  6  2  1  1  1  1  1  1 
summary(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2))
#   1  2  3  5  7  9 
# 17  8  5  1  3  1 
summary(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2))
#  1 2 
#  6 1

### PER PART ###
trip_count_median_1_16 <- favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_1_16$partid)$median
favstats(trip_count_median_1_16)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  1   2  1.2 0.4472 5       0
trip_count_median_17_34 <- favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_17_34$partid)$median
favstats(trip_count_median_17_34)
# min Q1 median Q3 max  mean  sd  n missing
#  1  1      1  2   4 1.619 0.9207 21       0
trip_count_median_35_54 <- favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_35_54$partid)$median
favstats(trip_count_median_35_54)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      2  3   7 2.294 1.611 17       0
trip_count_median_55_up <- favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_55_up$partid)$median
favstats(trip_count_median_55_up)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  1   1    1  0 4       0

trips_part2_age <- as.data.frame(cbind(c(trip_count_median_1_16, trip_count_median_17_34, trip_count_median_35_54, trip_count_median_55_up),
                                       c(rep("1-16",length(trip_count_median_1_16)),
                                         rep("17-34",length(trip_count_median_17_34)),
                                         rep("35-54",length(trip_count_median_35_54)),
                                         rep("55+",length(trip_count_median_55_up)))))
colnames(trips_part2_age) <- c("values","type")
trips_part2_age$values <- as.numeric(as.character(trips_part2_age$values))
trips_part2_age$type <- factor(trips_part2_age$type)

ggplot(data=trips_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.3) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,8.25)) + 
  scale_x_continuous(breaks=c(seq(0,10,2)))


## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,30))+ 
  scale_x_continuous(breaks=c(seq(0,30,2)))
#
##### percent biting time spent at home vs elsewhere #####
Nchelenge5_clusts_shorter2_biting_places3b_1_16 <- Nchelenge5_clusts_shorter2_biting_places3b[which(Nchelenge5_clusts_shorter2_biting_places3b$partid %in% Nchelenge_parts_1_16),]
Nchelenge5_clusts_shorter2_biting_places3b_17_34 <- Nchelenge5_clusts_shorter2_biting_places3b[which(Nchelenge5_clusts_shorter2_biting_places3b$partid %in% Nchelenge_parts_17_34),]
Nchelenge5_clusts_shorter2_biting_places3b_35_54 <- Nchelenge5_clusts_shorter2_biting_places3b[which(Nchelenge5_clusts_shorter2_biting_places3b$partid %in% Nchelenge_parts_35_54),]
Nchelenge5_clusts_shorter2_biting_places3b_55_up <- Nchelenge5_clusts_shorter2_biting_places3b[which(Nchelenge5_clusts_shorter2_biting_places3b$partid %in% Nchelenge_parts_55_up),]

favstats(Nchelenge5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time)
favstats(Nchelenge5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time[which(Nchelenge5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time != 0)])
favstats(Nchelenge5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time)
favstats(Nchelenge5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time[which(Nchelenge5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time != 0)])
favstats(Nchelenge5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time)
favstats(Nchelenge5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time[which(Nchelenge5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time != 0)])
favstats(Nchelenge5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time)
favstats(Nchelenge5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time[which(Nchelenge5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time != 0)])


perc_time_home2_age <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time, Nchelenge5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time,
                                             Nchelenge5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time, Nchelenge5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time),
                                           c(rep("1-16",length(Nchelenge5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time)),
                                             rep("17-34",length(Nchelenge5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time)),
                                             rep("35-54",length(Nchelenge5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time)),
                                             rep("55+",length(Nchelenge5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time)))))
colnames(perc_time_home2_age) <- c("values","type")
perc_time_home2_age$values <- as.numeric(as.character(perc_time_home2_age$values))
perc_time_home2_age$type <- factor(perc_time_home2_age$type)

ggplot(data=perc_time_home2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=10) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-1,101))
##
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
###### ###### RAINY/DRY SEASON SPLIT ###### ######
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
###############
###############
##### CHOMA #####
Choma5_new$rainy <- factor(0, levels=c(0,1))
Choma5_new$dry <- factor(0, levels=c(0,1))
Choma5_new$month <- as.numeric(c(unlist(strsplit(Choma5_new$date_new,"-"))[c(seq(from=2,to=c((nrow(Choma5_new)*3)-1), by=3))]))
Choma5_new$rainy[which(Choma5_new$month %in% c(1,2,3,10,11,12))] <- 1
Choma5_new$dry[which(Choma5_new$month %in% c(4:9))] <- 1
table(Choma5_new$dry,Choma5_new$partid) ## 15287 is only part with data in both rainy and dry
dry <- which(colSums(table(Choma5_new$dry,Choma5_new$partid)) == table(Choma5_new$dry,Choma5_new$partid)[2,]) ## 37 parts with data only in 1 season
rainy <- which(colSums(table(Choma5_new$rainy,Choma5_new$partid)) == table(Choma5_new$rainy,Choma5_new$partid)[2,]) ## 38 parts with data in both seasons
both <- which(!(1:nlevels(as.factor(Choma5_new$partid)) %in% c(dry,rainy))) ## one part with data in both seasons -- 15287

Choma5_clusts_shorter2$rainy <- factor(0, levels=c(0,1))
rainy <- Choma5_new$rainy[!duplicated(Choma5_new$partid)]
parts <- levels(as.factor(Choma5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Choma5_clusts_shorter2$rainy[which(Choma5_clusts_shorter2$partid == part)] <- rainy[ii]
}

Choma5_clusts_shorter2_home$rainy <- factor(0, levels=c(0,1))
rainy <- Choma5_new$rainy[!duplicated(Choma5_new$partid)]
parts <- levels(as.factor(Choma5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Choma5_clusts_shorter2_home$rainy[which(Choma5_clusts_shorter2_home$partid == part)] <- rainy[ii]
}

Choma5_clusts_shorter2_no_home$rainy <- factor(0, levels=c(0,1))
rainy <- Choma5_new$rainy[!duplicated(Choma5_new$partid)]
parts <- levels(as.factor(Choma5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Choma5_clusts_shorter2_no_home$rainy[which(Choma5_clusts_shorter2_no_home$partid == part)] <- rainy[ii]
}

Choma5_clust_trips_shorter2_no_home$rainy <- factor(0, levels=c(0,1))
rainy <- Choma5_new$rainy[!duplicated(Choma5_new$partid)]
parts <- levels(as.factor(Choma5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Choma5_clust_trips_shorter2_no_home$rainy[which(Choma5_clust_trips_shorter2_no_home$partid == part)] <- rainy[ii]
}

Choma5_clust_trips_shorter2_home$rainy <- factor(0, levels=c(0,1))
rainy <- Choma5_new$rainy[!duplicated(Choma5_new$partid)]
parts <- levels(as.factor(Choma5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Choma5_clust_trips_shorter2_home$rainy[which(Choma5_clust_trips_shorter2_home$partid == part)] <- rainy[ii]
}

Choma5_clust_trips_shorter2_home3$rainy <- factor(0, levels=c(0,1))
rainy <- Choma5_new$rainy[!duplicated(Choma5_new$partid)]
parts <- levels(as.factor(Choma5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Choma5_clust_trips_shorter2_home3$rainy[which(Choma5_clust_trips_shorter2_home3$partid == part)] <- rainy[ii]
}

Choma5_clust_trips_shorter2_no_home3$rainy <- factor(0, levels=c(0,1))
rainy <- Choma5_new$rainy[!duplicated(Choma5_new$partid)]
parts <- levels(as.factor(Choma5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Choma5_clust_trips_shorter2_no_home3$rainy[which(Choma5_clust_trips_shorter2_no_home3$partid == part)] <- rainy[ii]
}

Choma5_clust_trips_shorter2_home4$rainy <- factor(0, levels=c(0,1))
rainy <- Choma5_new$rainy[!duplicated(Choma5_new$partid)]
parts <- levels(as.factor(Choma5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Choma5_clust_trips_shorter2_home4$rainy[which(Choma5_clust_trips_shorter2_home4$partid == part)] <- rainy[ii]
}
#
#### FOR NOW, REMOVE 1 PART WITH RAINY AND DRY DATA ####
Choma5_rain_b <- Choma5_new[which(Choma5_new$partid != "15287"),]
Choma5_rain_clusts_shorter2 <- Choma5_clusts_shorter2[which(Choma5_clusts_shorter2$partid != "15287"),]
Choma5_rain_clusts_shorter2_home <- Choma5_clusts_shorter2_home[which(Choma5_clusts_shorter2_home$partid != "15287"),]
Choma5_rain_clusts_shorter2_no_home <- Choma5_clusts_shorter2_no_home[which(Choma5_clusts_shorter2_no_home$partid != "15287"),]
Choma5_rain_clust_trips_shorter2_no_home <- Choma5_clust_trips_shorter2_no_home[which(Choma5_clust_trips_shorter2_no_home$partid != "15287"),]
Choma5_rain_clust_trips_shorter2_home <- Choma5_clust_trips_shorter2_home[which(Choma5_clust_trips_shorter2_home$partid != "15287"),]
Choma5_rain_clust_trips_shorter2_home3 <- Choma5_clust_trips_shorter2_home3[which(Choma5_clust_trips_shorter2_home3$partid != "15287"),]
Choma5_rain_clust_trips_shorter2_no_home3 <- Choma5_clust_trips_shorter2_no_home3[which(Choma5_clust_trips_shorter2_no_home3$partid != "15287"),]
Choma5_rain_clust_trips_shorter2_home4 <- Choma5_clust_trips_shorter2_home4[which(Choma5_clust_trips_shorter2_home4$partid != "15287"),]
#
###### Split by Season ######
### number of locations (without home) ####
Choma5_rain_clusts_shorter2_no_home_dry <- Choma5_rain_clusts_shorter2_no_home[which(Choma5_rain_clusts_shorter2_no_home$rainy == 0),]
Choma5_rain_clusts_shorter2_no_home_rainy <- Choma5_rain_clusts_shorter2_no_home[which(Choma5_rain_clusts_shorter2_no_home$rainy == 1),]
loc_counts_dry <- favstats(Choma5_rain_clusts_shorter2_no_home_dry$clust~Choma5_rain_clusts_shorter2_no_home_dry$partid)$n
loc_counts_rainy <- favstats(Choma5_rain_clusts_shorter2_no_home_rainy$clust~Choma5_rain_clusts_shorter2_no_home_rainy$partid)$n
only_home_dry <- nlevels(as.factor(Choma5_rain_clusts_shorter2$partid[which(Choma5_rain_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Choma5_rain_clusts_shorter2_no_home_dry$partid))
only_home_rainy <- nlevels(as.factor(Choma5_rain_clusts_shorter2$partid[which(Choma5_rain_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Choma5_rain_clusts_shorter2_no_home_rainy$partid))
loc_counts_dry2 <- c(loc_counts_dry,rep(0, only_home_dry))
loc_counts_rainy2 <- c(loc_counts_rainy,rep(0, only_home_rainy))
favstats(loc_counts_dry2)
# min   Q1 median     Q3  max     mean       sd  n missing
#  0 6.75   10.5 17  20 11.21429 5.921443 28       0
favstats(loc_counts_rainy2)
# min Q1 median Q3 max     mean       sd  n missing
#   5  8     12 18  22 12.75758 5.420046 33       0

locs_count2 <- as.data.frame(cbind(c(loc_counts_dry2, loc_counts_rainy2),
                                   c(rep("dry",length(loc_counts_dry2)),
                                     rep("rainy",length(loc_counts_rainy2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

shapiro.test(locs_count2$values) ## not normal
kruskal.test(locs_count2$values, locs_count2$type) # p = 0.28

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-2,26)) +
  scale_x_continuous(breaks=c(seq(0,26,5)))
#

### distance of locations (without home) ####
### TOTAL ###
favstats(Choma5_rain_clusts_shorter2_no_home_dry$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#  0.03054886 0.6308518 1.635599 5.518482 112.2062 8.645072 18.65446 314       0
favstats(Choma5_rain_clusts_shorter2_no_home_rainy$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 0.07028766 6.764959 13.33677 17.63638 46.75496 12.41498 7.637024 421       0

### PER PART ###
medians_km_dry <- (favstats(Choma5_rain_clusts_shorter2_no_home_dry$hhdist_km_haversine~Choma5_rain_clusts_shorter2_no_home_dry$partid)$median)
medians_km_rainy <- (favstats(Choma5_rain_clusts_shorter2_no_home_rainy$hhdist_km_haversine~Choma5_rain_clusts_shorter2_no_home_rainy$partid)$median)
favstats(medians_km_dry) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
# 0.2448758 0.9144778 1.20238 4.978487 32.16992 4.399077 6.932997 27       0
favstats(medians_km_rainy) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.427719 1.368571 12.23068 16.42689 22.63018 10.9427 7.596466 33       0

## REMOVE OUTLIER ## none
medians_km_dry_short <- medians_km_dry
medians_km_rainy_short <- medians_km_rainy

dists_km2_season <- as.data.frame(cbind(c(medians_km_dry_short, medians_km_rainy_short),
                                        c(rep("dry",length(medians_km_dry_short)),
                                          rep("rainy",length(medians_km_rainy_short)))))
colnames(dists_km2_season) <- c("values","type")
dists_km2_season$values <- as.numeric(as.character(dists_km2_season$values))
dists_km2_season$type <- factor(dists_km2_season$type)

ggplot(data=dists_km2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,32))


kruskal.test(dists_km2_season$values, dists_km2_season$type)
# p = 0.002

#
#########
### distance of locations (without home) ####
Choma5_clusts_shorter2_no_home_1_16_rainy <- Choma5_clusts_shorter2_no_home_1_16[which(Choma5_clusts_shorter2_no_home_1_16$rainy == 1),]
Choma5_clusts_shorter2_no_home_17_34_rainy <- Choma5_clusts_shorter2_no_home_17_34[which(Choma5_clusts_shorter2_no_home_17_34$rainy == 1),]
Choma5_clusts_shorter2_no_home_35_54_rainy <- Choma5_clusts_shorter2_no_home_35_54[which(Choma5_clusts_shorter2_no_home_35_54$rainy == 1),]
Choma5_clusts_shorter2_no_home_55_up_rainy <- Choma5_clusts_shorter2_no_home_55_up[which(Choma5_clusts_shorter2_no_home_55_up$rainy == 1),]
Choma5_clusts_shorter2_no_home_1_16_dry <- Choma5_clusts_shorter2_no_home_1_16[which(Choma5_clusts_shorter2_no_home_1_16$rainy == 0),]
Choma5_clusts_shorter2_no_home_17_34_dry <- Choma5_clusts_shorter2_no_home_17_34[which(Choma5_clusts_shorter2_no_home_17_34$rainy == 0),]
Choma5_clusts_shorter2_no_home_35_54_dry <- Choma5_clusts_shorter2_no_home_35_54[which(Choma5_clusts_shorter2_no_home_35_54$rainy == 0),]
Choma5_clusts_shorter2_no_home_55_up_dry <- Choma5_clusts_shorter2_no_home_55_up[which(Choma5_clusts_shorter2_no_home_55_up$rainy == 0),]



medians_km_1_16_rainy <- (favstats(Choma5_clusts_shorter2_no_home_1_16_rainy$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_1_16_rainy$partid)$median)
medians_km_17_34_rainy <- (favstats(Choma5_clusts_shorter2_no_home_17_34_rainy$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_17_34_rainy$partid)$median)
medians_km_35_54_rainy <- (favstats(Choma5_clusts_shorter2_no_home_35_54_rainy$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_35_54_rainy$partid)$median)
medians_km_55_up_rainy <- (favstats(Choma5_clusts_shorter2_no_home_55_up_rainy$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_55_up_rainy$partid)$median)
medians_km_1_16_dry <- (favstats(Choma5_clusts_shorter2_no_home_1_16_dry$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_1_16_dry$partid)$median)
medians_km_17_34_dry <- (favstats(Choma5_clusts_shorter2_no_home_17_34_dry$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_17_34_dry$partid)$median)
medians_km_35_54_dry <- (favstats(Choma5_clusts_shorter2_no_home_35_54_dry$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_35_54_dry$partid)$median)
medians_km_55_up_dry <- (favstats(Choma5_clusts_shorter2_no_home_55_up_dry$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_55_up_dry$partid)$median)

dists_km2_age <- as.data.frame(cbind(c(medians_km_1_16_rainy, medians_km_17_34_rainy, medians_km_35_54_rainy, medians_km_55_up_rainy,
                                       medians_km_1_16_dry, medians_km_17_34_dry, medians_km_35_54_dry, medians_km_55_up_dry),
                                     c(rep("1-16",length(medians_km_1_16_rainy)),
                                       rep("17-34",length(medians_km_17_34_rainy)),
                                       rep("35-54",length(medians_km_35_54_rainy)),
                                       rep("55+",length(medians_km_55_up_rainy)),
                                       rep("1-16",length(medians_km_1_16_dry)),
                                       rep("17-34",length(medians_km_17_34_dry)),
                                       rep("35-54",length(medians_km_35_54_dry)),
                                       rep("55+",length(medians_km_55_up_dry))),
                                     c( rep("rainy", (length(medians_km_1_16_rainy) + length(medians_km_17_34_rainy) + length(medians_km_35_54_rainy) + length(medians_km_55_up_rainy))),
                                        rep("dry", (length(medians_km_1_16_dry) + length(medians_km_17_34_dry) + length(medians_km_35_54_dry) + length(medians_km_55_up_dry))))))


colnames(dists_km2_age) <- c("values","type","rainy")
dists_km2_age$values <- as.numeric(as.character(dists_km2_age$values))
dists_km2_age$type <- factor(dists_km2_age$type)
dists_km2_age$rainy <- factor(dists_km2_age$rainy)

favstats(dists_km2_age$values[which(dists_km2_age$rainy=="rainy")]~dists_km2_age$type[which(dists_km2_age$rainy=="rainy")])
favstats(dists_km2_age$values[which(dists_km2_age$rainy=="dry")]~dists_km2_age$type[which(dists_km2_age$rainy=="dry")])

dists_km2_age_rainy <- dists_km2_age[which(dists_km2_age$rainy == "rainy"),]
dists_km2_age_dry <- dists_km2_age[which(dists_km2_age$rainy == "dry"),]


kruskal.test(dists_km2_age_rainy$values, dists_km2_age_rainy$type) # p = 0.23
kruskal.test(dists_km2_age_dry$values, dists_km2_age_dry$type) # p = 0.23


pairwise.wilcox.test(dists_km2_age_rainy$values, dists_km2_age_rainy$type, p.adjust.method = "holm")
pairwise.wilcox.test(dists_km2_age_dry$values, dists_km2_age_rainy$type, p.adjust.method = "holm")

cc3 <- rbind(dists_km2_age_dry, dists_km2_age_rainy)
kruskal.test(cc3$values, cc3$type) # p = 0.23
pairwise.wilcox.test(cc3$values, cc3$type, p.adjust.method = "holm")


ggplot(data=dists_km2_age_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=12) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-3,65))


ggplot(data=dists_km2_age_dry, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=12) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-3,65))
###
#########
### number of trips (without home) ####
### TOTAL ###
Choma5_rain_clust_trips_shorter2_no_home3_dry <- Choma5_rain_clust_trips_shorter2_no_home3[which(Choma5_rain_clust_trips_shorter2_no_home3$rainy == 0),]
Choma5_rain_clust_trips_shorter2_no_home3_rainy <- Choma5_rain_clust_trips_shorter2_no_home3[which(Choma5_rain_clust_trips_shorter2_no_home3$rainy == 1),]

favstats(Choma5_rain_clust_trips_shorter2_no_home3_dry$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      2  3  52 3.009554 4.777801 314       0
favstats(Choma5_rain_clust_trips_shorter2_no_home3_rainy$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      1  3  58 3.106888 5.645097 421       0

### PER PART ###
trips_part_median_dry <- favstats(Choma5_rain_clust_trips_shorter2_no_home3_dry$trip_count_new ~ Choma5_rain_clust_trips_shorter2_no_home3_dry$partid)$median
trips_part_median_rainy <- favstats(Choma5_rain_clust_trips_shorter2_no_home3_rainy$trip_count_new ~ Choma5_rain_clust_trips_shorter2_no_home3_rainy$partid)$median
favstats(trips_part_median_dry)
#    min Q1 median  Q3   max      mean       sd  n missing
#  1    1      1  2   4 1.592593 0.7970744 27       0
favstats(trips_part_median_rainy)
#   min Q1 median  Q3   max      mean       sd  n missing
#     1  1    1.5  2   3 1.515152 0.5517706 33       0

trips_part2_season <- as.data.frame(cbind(c(trips_part_median_dry, trips_part_median_rainy),
                                          c(rep("dry",length(trips_part_median_dry)),
                                            rep("rainy",length(trips_part_median_rainy)))))
colnames(trips_part2_season) <- c("values","type")
trips_part2_season$values <- as.numeric(as.character(trips_part2_season$values))
trips_part2_season$type <- factor(trips_part2_season$type)

ggplot(data=trips_part2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks=c(0:10))
#
kruskal.test(trips_part2_season$values, trips_part2_season$type)
# p = 0.97

### time per trip ####
### TOTAL ###
Choma5_rain_clust_trips_shorter2_no_home_dry <- Choma5_rain_clust_trips_shorter2_no_home[which(Choma5_rain_clust_trips_shorter2_no_home$rainy == 0),]
Choma5_rain_clust_trips_shorter2_no_home_rainy <- Choma5_rain_clust_trips_shorter2_no_home[which(Choma5_rain_clust_trips_shorter2_no_home$rainy == 1),]
favstats(Choma5_rain_clust_trips_shorter2_no_home_dry$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2402778 0.6183333 1.2525 3.095556 232.08 4.468468 13.35561 931       0
favstats(Choma5_rain_clust_trips_shorter2_no_home_rainy$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2402778 0.7197222 1.554861 5.230833 472.56 11.85264 38.64956 1298       0

### PER PART ###
trip_time_part_median_dry <- favstats(Choma5_rain_clust_trips_shorter2_no_home_dry$trip_time_new_hour ~ Choma5_rain_clust_trips_shorter2_no_home_dry$partid)$median
trip_time_part_median_rainy <- favstats(Choma5_rain_clust_trips_shorter2_no_home_rainy$trip_time_new_hour ~ Choma5_rain_clust_trips_shorter2_no_home_rainy$partid)$median
favstats(trip_time_part_median_dry)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.5879167 0.9297222 1.286111 1.526181 3.36 1.378524 0.6197525 27       0
favstats(trip_time_part_median_rainy)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.6075 1.221389 1.474444 1.752222 5.650556 1.773506 1.073157 33       0

## REMOVE OUTLIERS ##
trip_time_part_median_dry_short <- trip_time_part_median_dry ## no outliers to remove

trip_times_part2_season <- as.data.frame(cbind(c(trip_time_part_median_dry_short, trip_time_part_median_rainy),
                                               c(rep("dry",length(trip_time_part_median_dry_short)),
                                                 rep("rainy",length(trip_time_part_median_rainy)))))
colnames(trip_times_part2_season) <- c("values","type")
trip_times_part2_season$values <- as.numeric(as.character(trip_times_part2_season$values))
trip_times_part2_season$type <- factor(trip_times_part2_season$type)

shapiro.test(trip_times_part2_season$values) ## not normal
kruskal.test(trip_times_part2_season$values, trip_times_part2_season$type) # p = 0.08

ggplot(data=trip_times_part2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 0.5, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,6)) +
  scale_x_continuous(breaks=c(0:6))
#
### Home Loc ####
Choma5_rain_clusts_shorter2_home_dry <- Choma5_rain_clusts_shorter2_home[which(Choma5_rain_clusts_shorter2_home$rainy == 0),] 
Choma5_rain_clusts_shorter2_home_rainy <- Choma5_rain_clusts_shorter2_home[which(Choma5_rain_clusts_shorter2_home$rainy == 1),] 

no_home_dry <- nlevels(as.factor(Choma5_rain_clusts_shorter2$partid[which(Choma5_rain_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Choma5_rain_clusts_shorter2_home_dry$partid))
## 3 drys with no home
no_home_rainy <- nlevels(as.factor(Choma5_rain_clusts_shorter2$partid[which(Choma5_rain_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Choma5_rain_clusts_shorter2_home_rainy$partid))
## no rainys with no home

### number of 'trips' to home loc ####
### TOTAL ###
Choma5_rain_clust_trips_shorter2_home3_dry <- Choma5_rain_clust_trips_shorter2_home3[which(Choma5_rain_clust_trips_shorter2_home3$rainy == 0),]
Choma5_rain_clust_trips_shorter2_home3_rainy <- Choma5_rain_clust_trips_shorter2_home3[which(Choma5_rain_clust_trips_shorter2_home3$rainy == 1),]

favstats(Choma5_rain_clust_trips_shorter2_home3_dry$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#  2 5.75   11.5 30.25  38 16.41667 12.87369 24       0
favstats(Choma5_rain_clust_trips_shorter2_home3_rainy$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#   1  4      4 11  35 7.757576 7.314499 33       0

### total time at home loc ####
Choma5_rain_clust_trips_shorter2_home4_dry <- Choma5_rain_clust_trips_shorter2_home4[which(Choma5_rain_clust_trips_shorter2_home4$rainy == 0),]
Choma5_rain_clust_trips_shorter2_home4_rainy <- Choma5_rain_clust_trips_shorter2_home4[which(Choma5_rain_clust_trips_shorter2_home4$rainy == 1),]

favstats(Choma5_rain_clust_trips_shorter2_home4_dry$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
# 0.5471991 20.11667 24.12972 27.22779 32.9477 22.82372 7.589278 24       0
favstats(Choma5_rain_clust_trips_shorter2_home4_rainy$clust_trips_time_new)
#         min       Q1   median       Q3      max     mean       sd  n missing
# 0.02523149 2.416597 7.155851 25.68119 40.365 14.37288 13.225 33       0
#
### percent time at home ####
favstats(Choma5_rain_clusts_shorter2_home_dry$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
#  2.184087 82.40033 89.01833 97.44476 100 83.88992 23.0591 25       0
favstats(Choma5_rain_clusts_shorter2_home_rainy$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#  0.08520987 7.064713 19.4021 90.47754 98.82904 41.19878 39.22825 33       0

perc_time_home2_season <- as.data.frame(cbind(c(Choma5_rain_clusts_shorter2_home_dry$perc_time_home, Choma5_rain_clusts_shorter2_home_rainy$perc_time_home),
                                              c(rep("dry",length(Choma5_rain_clusts_shorter2_home_dry$perc_time_home)),
                                                rep("rainy",length(Choma5_rain_clusts_shorter2_home_rainy$perc_time_home)))))
colnames(perc_time_home2_season) <- c("values","type")
perc_time_home2_season$values <- as.numeric(as.character(perc_time_home2_season$values))
perc_time_home2_season$type <- factor(perc_time_home2_season$type)

ggplot(data=perc_time_home2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.25) +
  scale_fill_manual("Season",values = c("sienna", "skyblue"), labels=c("Dry","Rainy")) +
  scale_color_manual("Season",values = c("sienna", "skyblue"), labels=c("Dry","Rainy")) +
  scale_x_continuous("Percent time at home", breaks=c(seq(0,100,25))) +
  coord_cartesian(xlim=c(0,100)) +
  theme_classic() +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=12),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12))

kruskal.test(perc_time_home2_season$values, perc_time_home2_season$type)
# p < 0.001

##
### Average Direction of Travel #####
Choma5_rain_clusts_shorter2_no_home$rainy <- as.factor(Choma5_rain_clusts_shorter2_no_home$rainy)
temp <- Choma5_rain_clusts_shorter2_no_home[!duplicated(Choma5_rain_clusts_shorter2_no_home$partid),]
temp2 <- temp[which(!(is.na(temp$rainy))),]

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=rainy), position=position_dodge2(preserve = "single"), width=3) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Rainy Season", values=c("sienna","skyblue")) +
  theme(legend.position = "bottom")
#

favstats(temp2$avg_dir_from_home_bearing~temp2$rainy)
# temp2$rainy      min       Q1    median        Q3      max     mean       sd  n missing
# 1           0 43.53911 123.6301 183.61081 279.8227 334.6366 200.7264 89.98696 27       0
# 2           1  8.43748  53.8835  73.97805 103.6810 344.5980 102.8886 89.70920 33       0
favstats(temp2$avg_dir_from_home_mag~temp2$rainy)
# temp2$rainy       min        Q1    median        Q3      max      mean       sd  n missing
# 1           0 0.1915475 0.7464408  1.649824  9.713094 33.03357  6.790245 9.648910 27       0
# 2           1 0.3610027 1.6651872 11.900858 14.832664 22.15730 10.167155 6.914351 33       0


temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)] <- -(360-temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)])
ggplot(data=temp2, aes(x=avg_dir_from_home_values, y=..count.., fill=rainy)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~rainy, labeller = labeller(rainy = c("0"="Dry", "1"="Rainy")))+
  scale_fill_manual(values=c("sienna","skyblue")) +
  scale_y_continuous(breaks=c(seq(0,10,2)))





####
temp2$perc_time_home <- as.numeric(NA)
temp2$perc_time_home[which(temp2$partid %in% levels(as.factor(Choma5_rain_clusts_shorter2_home_dry$partid)))] <- Choma5_rain_clusts_shorter2_home_dry$perc_time_home[which(Choma5_rain_clusts_shorter2_home_dry$partid %in% levels(as.factor(temp2$partid)))]
temp2$perc_time_home[which(temp2$partid %in% levels(as.factor(Choma5_rain_clusts_shorter2_home_rainy$partid)))] <- Choma5_rain_clusts_shorter2_home_rainy$perc_time_home[which(Choma5_rain_clusts_shorter2_home_rainy$partid %in% levels(as.factor(temp2$partid)))]
# 15287 --> rainy and dry --> not included for split
# 15128 --> only home --> no direction/dist
# "10692b" "10753"  "15927"  --> no home loc --> no percentage
temp3 <- temp2[which(!(is.na(temp2$perc_time_home))),]

ggplot(data=temp3, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=rainy, alpha=perc_time_home), position=position_dodge2(preserve = "single"), width=3) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Rainy Season", values=c("sienna","skyblue")) +
  theme(legend.position = "bottom")
#


## average direction calc
temp2$values2 <- 90 - temp2$avg_dir_from_home_values
temp2$values2[which(temp2$values2 < 0)] <- temp2$values2[which(temp2$values2 < 0)] + 360
temp2$values3 <- deg2rad(temp2$values2)

temp2$avg_value <- as.numeric(NA)
for(ii in 1:nlevels(temp2$rainy)){
  level <- levels(temp2$rainy)[ii]
  temp <- temp2[which(temp2$rainy == level),]
  sum_cos_dist <- sum(cos(temp$values3))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$values3))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  temp2$avg_value[which(temp2$rainy == level)] <- avg_deg
}

temp2$avg_value2 <- 90 - temp2$avg_value
temp2$avg_value2[which(temp2$avg_value2 < 0)] <- temp2$avg_value2[which(temp2$avg_value2 < 0)] + 360
favstats(temp2$avg_value2~temp2$rainy)
## dry: 187.14 --> S
## rainy: 53.6 --> NE

kruskal.test(temp2$avg_value2, temp2$rainy)


aa_rainy <- Choma5_rain_clusts_shorter2_no_home[which(Choma5_rain_clusts_shorter2_no_home$rainy == 1),]
aa_rainy2 <- aa_rainy[,c(1,8:11,15,19,40:50)]
bb <- aa_rainy2[!duplicated(aa_rainy2$partid),]
bb$home_clust2 <- 1
bb$clust <- 99
bb$log_lat <- bb$hh_lat
bb$log_long <- bb$hh_long
aa_rainy2b <- rbind(aa_rainy2,bb)
aa_rainy3 <- aa_rainy2b[order(aa_rainy2b$partid,-aa_rainy2b$home_clust2),]
aa_rainy3$home_west <- factor(0,levels=c(0,1))
aa_rainy3$home_west[which(aa_rainy3$hh_lat < -16.33)] <- 1
bb <- aa_rainy3[which(aa_rainy3$home_west == 1),]
aa2 <- aa_rainy3[which(aa_rainy3$partid == "10334"),]

xy <- SpatialPointsDataFrame(
  matrix(c(aa_rainy3$log_lat, aa_rainy3$log_long), ncol=2), data.frame(ID=seq(1:nrow(aa_rainy3))),
  proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
plot(xy, axes=TRUE, col=c(rainbow(33))[as.factor(aa_rainy3$partid)],
     pch=ifelse(aa_rainy3$home_clust2 == 0, 1, 8),
     cex=ifelse(aa_rainy3$home_clust2 == 1, 3.5, 0.75),
     xlim=c(-16.5,-16.2),ylim=c(26.78,26.98))


xy <- SpatialPointsDataFrame(
  matrix(c(aa_rainy3$log_lat, aa_rainy3$log_long), ncol=2), data.frame(ID=seq(1:nrow(aa_rainy3))),
  proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
plot(xy, axes=TRUE, col=c(rainbow(11),rep("black",2),rainbow(2),rep("black",2),rainbow(15),"black")[as.factor(aa_rainy3$partid)],
     pch=ifelse(aa_rainy3$home_clust2 == 0, 1, 8),
     cex=ifelse(aa_rainy3$home_clust2 == 1, 3.5, 0.75),
     xlim=c(-16.5,-16.2),ylim=c(26.78,26.98))


xy <- SpatialPointsDataFrame(
  matrix(c(bb$log_lat, bb$log_long), ncol=2), data.frame(ID=seq(1:nrow(bb))),
  proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
plot(xy, axes=TRUE, col=c(rainbow(17),rep("gray",3),rainbow(3),"gray",rainbow(1),"grey",rainbow(1),"grey")[as.factor(bb$partid)],
     pch=ifelse(bb$home_clust2 == 0, 1, 8),
     cex=ifelse(bb$home_clust2 == 1, 3.5, 0.75),
     xlim=c(-16.5,-16.2),ylim=c(26.78,26.98))

xy <- SpatialPointsDataFrame(
  matrix(c(bb$log_lat, bb$log_long), ncol=2), data.frame(ID=seq(1:nrow(bb))),
  proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
plot(xy, axes=TRUE, col=c(rep("gray",17),rainbow(1))[as.factor(bb$partid)],
     pch=ifelse(bb$home_clust2 == 0, 1, 8),
     cex=ifelse(bb$home_clust2 == 1, 3.5, 0.75),
     xlim=c(-16.5,-16.2),ylim=c(26.78,26.98))








xy <- SpatialPointsDataFrame(
  matrix(c(aa_rainy3$log_lat, aa_rainy3$log_long), ncol=2), data.frame(ID=seq(1:nrow(aa_rainy3))),
  proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
plot(xy, axes=TRUE, col=c(rainbow(11),rep("black",2),rainbow(2),rep("black",2),rainbow(15),"black")[as.factor(aa_rainy3$partid)],
     pch=ifelse(aa_rainy3$home_clust2 == 0, 1, 8),
     cex=ifelse(aa_rainy3$home_clust2 == 1, 3.5, 0.75),
     xlim=c(-16.5,-16.2),ylim=c(26.78,26.98))



sbbox <- make_bbox(lon = c(26.73,27.02), lat = c(-16.5,-16.23), f = .1)
# get map
choma = get_map(location=sbbox, maptype="terrain-labels",zoom=13)
# create map
chomamap2 = ggmap(choma)
# display map
chomamap2 +
  geom_point(data = aa_rainy3, mapping = aes(x = hh_long, y = hh_lat), 
             color = c(rainbow(12),rep("black",1),rainbow(2),rep("black",2),rainbow(1),
                       "black",rainbow(8),rep("black",3),rainbow(1),
                       "black",rep("black",2))[as.factor(aa_rainy3$partid)],
             shape=8,size=7)  +
  
geom_point(data = aa_rainy3, mapping = aes(x = log_long, y = log_lat), 
             color = c(rainbow(12),rep("black",1),rainbow(2),rep("black",2),rainbow(1),
                       "black",rainbow(8),rep("black",3),rainbow(1),
                       "black",rep("black",2))[as.factor(aa_rainy3$partid)],
           shape=1) 
  
  
#
### time spent at clusters vs elsewhere (travel, less important places, etc) ####
# percent time in clusters
Choma5_rain_clusts_shorter2_dry <- Choma5_rain_clusts_shorter2[which(Choma5_rain_clusts_shorter2$rainy == 0),]
Choma5_rain_clusts_shorter2_rainy <- Choma5_rain_clusts_shorter2[which(Choma5_rain_clusts_shorter2$rainy == 1),]
Choma5_rain_clusts_shorter2_dry$perc_time_in_all_clusts <- Choma5_rain_clusts_shorter2_dry$part_trips_time_new/Choma5_rain_clusts_shorter2_dry$total_time
Choma5_rain_clusts_shorter2_rainy$perc_time_in_all_clusts <- Choma5_rain_clusts_shorter2_rainy$part_trips_time_new/Choma5_rain_clusts_shorter2_rainy$total_time

summary(Choma5_rain_clusts_shorter2_dry$perc_time_in_all_clusts[!duplicated(Choma5_rain_clusts_shorter2_dry$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6629  0.8500  0.9261  0.8949  0.9534  0.9955 
summary(Choma5_rain_clusts_shorter2_rainy$perc_time_in_all_clusts[!duplicated(Choma5_rain_clusts_shorter2_rainy$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8287  0.8902  0.9373  0.9318  0.9765  0.9950 

# time outside clusters
summary(Choma5_rain_clusts_shorter2_dry$total_time[!duplicated(Choma5_rain_clusts_shorter2_dry$partid)] - Choma5_rain_clusts_shorter2_dry$part_trips_time_new[!duplicated(Choma5_rain_clusts_shorter2_dry$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1395  1.1869  2.2344  3.2449  4.6457  9.8070 
summary(Choma5_rain_clusts_shorter2_rainy$total_time[!duplicated(Choma5_rain_clusts_shorter2_rainy$partid)] - Choma5_rain_clusts_shorter2_rainy$part_trips_time_new[!duplicated(Choma5_rain_clusts_shorter2_rainy$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2383  0.8182  1.8176  2.4824  3.9561  6.6361 
#
#
#
### Time at locations ####
favstats(Choma5_rain_clusts_shorter2_no_home_dry$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2955556 1.091319 2.319722 4.821458 684.7197 13.26256 57.84376 314       0
favstats(Choma5_rain_clusts_shorter2_no_home_rainy$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2461111 1.0925 1.690556 4.628889 885.6553 36.55088 139.2314 421       0

### PER PART ###
clust_time_part_median_dry <- (favstats(Choma5_rain_clusts_shorter2_no_home_dry$clust_trips_time_new~Choma5_rain_clusts_shorter2_no_home_dry$partid)$median)*24
clust_time_part_median_rainy <- (favstats(Choma5_rain_clusts_shorter2_no_home_rainy$clust_trips_time_new~Choma5_rain_clusts_shorter2_no_home_rainy$partid)$median)*24
favstats(clust_time_part_median_dry)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.7801389 1.621389 1.986667 2.700694 6.960278 2.396404 1.348462 27       0
favstats(clust_time_part_median_rainy)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.08 1.369722   1.68 2.048056 7.662778 1.992113 1.171242 33       0

clust_time_part_rainy <- as.data.frame(cbind(c(clust_time_part_median_dry, clust_time_part_median_rainy),
                                              c(rep("dry",length(clust_time_part_median_dry)),
                                                rep("rainy",length(clust_time_part_median_rainy)))))
colnames(clust_time_part_rainy) <- c("values","type")
clust_time_part_rainy$values <- as.numeric(as.character(clust_time_part_rainy$values))
clust_time_part_rainy$type <- factor(clust_time_part_rainy$type)

shapiro.test(clust_time_part_rainy$values) ## not normal
kruskal.test(clust_time_part_rainy$values, clust_time_part_rainy$type) # p = 0.12

ggplot(data=clust_time_part_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,9)) +
  scale_x_continuous(breaks=c(seq(0,9,1)))
#
#
#
####### Biting time metrics #####
### number locations outside home ####
Choma5_clusts_shorter2_no_home_biting_places <- Choma5_clusts_shorter2_no_home[which(!(is.na(Choma5_clusts_shorter2_no_home$clust_trips_time_new_biting_time))),]
Choma5_clusts_shorter2_no_home_biting_places2 <- Choma5_clusts_shorter2_no_home_biting_places[which((Choma5_clusts_shorter2_no_home_biting_places$clust_trips_time_new_biting_time*24) >= 0.5),]

Choma5_clusts_shorter2_no_home_biting_places2_dry <- Choma5_clusts_shorter2_no_home_biting_places2[which(Choma5_clusts_shorter2_no_home_biting_places2$rainy == 0),]
Choma5_clusts_shorter2_no_home_biting_places2_rainy <- Choma5_clusts_shorter2_no_home_biting_places2[which(Choma5_clusts_shorter2_no_home_biting_places2$rainy == 1),]

loc_counts_dry <- favstats(Choma5_clusts_shorter2_no_home_biting_places2_dry$clust~Choma5_clusts_shorter2_no_home_biting_places2_dry$partid)$n
loc_counts_rainy <- favstats(Choma5_clusts_shorter2_no_home_biting_places2_rainy$clust~Choma5_clusts_shorter2_no_home_biting_places2_rainy$partid)$n
only_home_dry <- nlevels(as.factor(Choma5_rain_clusts_shorter2$partid[which(Choma5_rain_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_dry$partid))
only_home_rainy <- nlevels(as.factor(Choma5_rain_clusts_shorter2$partid[which(Choma5_rain_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_rainy$partid))
loc_counts_dry2 <- c(loc_counts_dry,rep(0, only_home_dry))
loc_counts_rainy2 <- c(loc_counts_rainy,rep(0, only_home_rainy))

favstats(loc_counts_dry2)
# min Q1 median Q3 max     mean       sd  n missing
#  0  0      1  3   5 1.571 1.665 28       0
favstats(loc_counts_rainy2)
# min Q1 median Q3 max     mean       sd  n missing
#    0  1      2  3   6 2.273 1.825 33       0
summary(as.factor(loc_counts_dry2))
#  0  1  2  3  4  5 
# 11  5  3  5  2  2 
summary(as.factor(loc_counts_rainy2))
# 0 1 2 3 4 5 6 
# 5 9 6 6 2 2 3 
## all those with 0 have biting time at home location (not just nowhere)

locs_count2 <- as.data.frame(cbind(c(loc_counts_dry2, loc_counts_rainy2),
                                   c(rep("dry",length(loc_counts_dry2)),
                                     rep("rainy",length(loc_counts_rainy2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) + 
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-0.35,6.35)) +
  scale_x_continuous(breaks=c(seq(0,7,1)))
#
##
### time per location ####
favstats((Choma5_clusts_shorter2_no_home_biting_places2_dry$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5383 1.998  4.344 18.89 218.4 25.94 53.01 44       0
favstats((Choma5_clusts_shorter2_no_home_biting_places2_rainy$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5142 1.237  6.614 137.5 351.6 62.77 94.79 75       0

clusts_median_Choma_dry <- (favstats(Choma5_clusts_shorter2_no_home_biting_places2_dry$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_dry$partid)$median)*24
clusts_median_Choma_rainy <- (favstats(Choma5_clusts_shorter2_no_home_biting_places2_rainy$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_rainy$partid)$median)*24
favstats(clusts_median_Choma_dry)
# min       Q1   median  Q3      max     mean       sd  n missing
# 1.67 2.16  3.968 17.4 218.4 35.15 65.81 17       0
favstats(clusts_median_Choma_rainy)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.5928 3.878  12.63 105.2 291.3 64.02 81.91 28       0

clust_time_part_rainy <- as.data.frame(cbind(c(clusts_median_Choma_dry, clusts_median_Choma_rainy),
                                              c(rep("dry",length(clusts_median_Choma_dry)),
                                                rep("rainy",length(clusts_median_Choma_rainy)))))
colnames(clust_time_part_rainy) <- c("values","type")
clust_time_part_rainy$values <- as.numeric(as.character(clust_time_part_rainy$values))
clust_time_part_rainy$type <- factor(clust_time_part_rainy$type)

ggplot(data=clust_time_part_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=2.5) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,300)) +
  scale_x_continuous(breaks=c(seq(0,300,50)))
#
### number trips ####
Choma5_clust_trips_shorter2_no_home_biting_places <- Choma5_clust_trips_shorter2_no_home[which(!(is.na(Choma5_clust_trips_shorter2_no_home$trip_time_new2_biting_time))),]
## only when >30 minutes during biting time
Choma5_clust_trips_shorter2_no_home_biting_places2 <- Choma5_clust_trips_shorter2_no_home_biting_places[which((Choma5_clust_trips_shorter2_no_home_biting_places$trip_time_new2_biting_time*24) >= 0.5),]
## removes 63 trips, 3 participants

Choma5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2 <- numeric(length=nrow(Choma5_clust_trips_shorter2_no_home_biting_places2))
parts <- levels(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places2$clust[which(Choma5_clust_trips_shorter2_no_home_biting_places2$partid == part)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    counts <- Choma5_clust_trips_shorter2_no_home_biting_places2$trip_count_new[which(Choma5_clust_trips_shorter2_no_home_biting_places2$partid == part & Choma5_clust_trips_shorter2_no_home_biting_places2$clust == clust)]
    counts <- as.factor(counts)
    levels(counts) <- c(1:nlevels(counts))
    Choma5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2[which(Choma5_clust_trips_shorter2_no_home_biting_places2$partid == part & Choma5_clust_trips_shorter2_no_home_biting_places2$clust == clust)] <- as.numeric(as.character(counts))
  }
}

Choma5_clust_trips_shorter2_no_home_biting_places2b <- Choma5_clust_trips_shorter2_no_home_biting_places2[order(Choma5_clust_trips_shorter2_no_home_biting_places2$partid,
                                                                                                                Choma5_clust_trips_shorter2_no_home_biting_places2$clust,
                                                                                                                -Choma5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2),]
Choma5_clust_trips_shorter2_no_home_biting_places3 <- Choma5_clust_trips_shorter2_no_home_biting_places2b[!duplicated(Choma5_clust_trips_shorter2_no_home_biting_places2b[,c(1,15)]),]

Choma5_clust_trips_shorter2_no_home_biting_places3_dry <- Choma5_clust_trips_shorter2_no_home_biting_places3[which(Choma5_clust_trips_shorter2_no_home_biting_places3$rainy == 0),]
Choma5_clust_trips_shorter2_no_home_biting_places3_rainy <- Choma5_clust_trips_shorter2_no_home_biting_places3[which(Choma5_clust_trips_shorter2_no_home_biting_places3$rainy == 1),]

favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2  27 3.326 5.931 43       0
favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#  1  1      2 4.5  27 4.493 5.999 75       0

summary(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2))
#   1  2  3  4  5  8 22 23 27 
#  26  7  3  1  2  1  1  1  1 
summary(as.factor(Choma5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2))
#  1  2  3  4  5  6  7  8 10 12 13 14 15 17 21 22 27 
# 36 13  5  2  1  1  1  2  3  2  2  1  1  1  1  2  1 

### PER PART ###
trip_count_median_dry <- favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_dry$partid)$median
favstats(trip_count_median_dry)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1    1.5  2  27 4.353 7.691 17       0
trip_count_median_rainy <- favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_rainy$partid)$median
favstats(trip_count_median_rainy)
# min Q1 median Q3 max  mean  sd  n missing
#   1   1      2 4.5  12 3.143 2.748 28       0

trips_part2_rainy <- as.data.frame(cbind(c(trip_count_median_dry, trip_count_median_rainy),
                                          c(rep("dry",length(trip_count_median_dry)),
                                            rep("rainy",length(trip_count_median_rainy)))))
colnames(trips_part2_rainy) <- c("values","type")
trips_part2_rainy$values <- as.numeric(as.character(trips_part2_rainy$values))
trips_part2_rainy$type <- factor(trips_part2_rainy$type)

ggplot(data=trips_part2_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,26.25)) + 
  scale_x_continuous(breaks=c(seq(0,30,2)))


## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,30))+ 
  scale_x_continuous(breaks=c(seq(0,30,2)))
#
##### percent biting time spent at home vs elsewhere #####
Choma5_clusts_shorter2_biting_places <- Choma5_clusts_shorter2[which(!(is.na(Choma5_clusts_shorter2$clust_trips_time_new_biting_time))),]
Choma5_clusts_shorter2_biting_places$percent_home_biting_time <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_biting_places$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  percent <- Choma5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Choma5_clusts_shorter2_biting_places$partid == part & Choma5_clusts_shorter2_biting_places$home_clust2 == 1)]/
    sum(Choma5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Choma5_clusts_shorter2_biting_places$partid == part)])
  if(length(percent) == 0){
    percent <- 0
  }
  Choma5_clusts_shorter2_biting_places$percent_home_biting_time[which(Choma5_clusts_shorter2_biting_places$partid == part)] <- percent*100
}
Choma5_clusts_shorter2_biting_places2 <- Choma5_clusts_shorter2_biting_places[!duplicated(Choma5_clusts_shorter2_biting_places$partid),]
homes <- Choma5_clusts_shorter2_home$partid
Choma5_clusts_shorter2_biting_places3 <- Choma5_clusts_shorter2_biting_places[which(Choma5_clusts_shorter2_biting_places$partid %in% homes),]
Choma5_clusts_shorter2_biting_places3b <- Choma5_clusts_shorter2_biting_places3[!duplicated(Choma5_clusts_shorter2_biting_places3$partid),]

Choma5_clusts_shorter2_biting_places3b_dry <- Choma5_clusts_shorter2_biting_places3b[which(Choma5_clusts_shorter2_biting_places3b$rainy == 0),]
Choma5_clusts_shorter2_biting_places3b_rainy <- Choma5_clusts_shorter2_biting_places3b[which(Choma5_clusts_shorter2_biting_places3b$rainy == 1),]

favstats(Choma5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time)
favstats(Choma5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time[which(Choma5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time != 0)])
favstats(Choma5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time)
favstats(Choma5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time[which(Choma5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time != 0)])


perc_time_home2_rainy <- as.data.frame(cbind(c(Choma5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time, Choma5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time),
                                              c(rep("dry",length(Choma5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time)),
                                                rep("rainy",length(Choma5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time)))))
colnames(perc_time_home2_rainy) <- c("values","type")
perc_time_home2_rainy$values <- as.numeric(as.character(perc_time_home2_rainy$values))
perc_time_home2_rainy$type <- factor(perc_time_home2_rainy$type)

ggplot(data=perc_time_home2_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-1,101))
##
###############
###############
##### MUTASA #####
Mutasa5_new$rainy <- factor(0, levels=c(0,1))
Mutasa5_new$dry <- factor(0, levels=c(0,1))
Mutasa5_new$month <- as.numeric(c(unlist(strsplit(Mutasa5_new$date_new,"-"))[c(seq(from=2,to=c((nrow(Mutasa5_new)*3)-1), by=3))]))
Mutasa5_new$rainy[which(Mutasa5_new$month %in% c(1,2,3,10,11,12))] <- 1
Mutasa5_new$dry[which(Mutasa5_new$month %in% c(4:9))] <- 1
table(Mutasa5_new$dry,Mutasa5_new$partid) ## 15287 is only part with data in both rainy and dry
dry <- which(colSums(table(Mutasa5_new$dry,Mutasa5_new$partid)) == table(Mutasa5_new$dry,Mutasa5_new$partid)[2,]) ## 37 parts with data only in 1 season
rainy <- which(colSums(table(Mutasa5_new$rainy,Mutasa5_new$partid)) == table(Mutasa5_new$rainy,Mutasa5_new$partid)[2,]) ## 38 parts with data in both seasons
both <- which(!(1:nlevels(as.factor(Mutasa5_new$partid)) %in% c(dry,rainy))) ## no parts with data in both seasons

Mutasa5_clusts_shorter2$rainy <- factor(0, levels=c(0,1))
rainy <- Mutasa5_new$rainy[!duplicated(Mutasa5_new$partid)]
parts <- levels(as.factor(Mutasa5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Mutasa5_clusts_shorter2$rainy[which(Mutasa5_clusts_shorter2$partid == part)] <- rainy[ii]
}

Mutasa5_clusts_shorter2_home$rainy <- factor(0, levels=c(0,1))
rainy <- Mutasa5_new$rainy[!duplicated(Mutasa5_new$partid)]
parts <- levels(as.factor(Mutasa5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Mutasa5_clusts_shorter2_home$rainy[which(Mutasa5_clusts_shorter2_home$partid == part)] <- rainy[ii]
}

Mutasa5_clusts_shorter2_no_home$rainy <- factor(0, levels=c(0,1))
rainy <- Mutasa5_new$rainy[!duplicated(Mutasa5_new$partid)]
parts <- levels(as.factor(Mutasa5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Mutasa5_clusts_shorter2_no_home$rainy[which(Mutasa5_clusts_shorter2_no_home$partid == part)] <- rainy[ii]
}

Mutasa5_clust_trips_shorter2_no_home$rainy <- factor(0, levels=c(0,1))
rainy <- Mutasa5_new$rainy[!duplicated(Mutasa5_new$partid)]
parts <- levels(as.factor(Mutasa5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Mutasa5_clust_trips_shorter2_no_home$rainy[which(Mutasa5_clust_trips_shorter2_no_home$partid == part)] <- rainy[ii]
}

Mutasa5_clust_trips_shorter2_home$rainy <- factor(0, levels=c(0,1))
rainy <- Mutasa5_new$rainy[!duplicated(Mutasa5_new$partid)]
parts <- levels(as.factor(Mutasa5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Mutasa5_clust_trips_shorter2_home$rainy[which(Mutasa5_clust_trips_shorter2_home$partid == part)] <- rainy[ii]
}

Mutasa5_clust_trips_shorter2_home3$rainy <- factor(0, levels=c(0,1))
rainy <- Mutasa5_new$rainy[!duplicated(Mutasa5_new$partid)]
parts <- levels(as.factor(Mutasa5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Mutasa5_clust_trips_shorter2_home3$rainy[which(Mutasa5_clust_trips_shorter2_home3$partid == part)] <- rainy[ii]
}

Mutasa5_clust_trips_shorter2_no_home3$rainy <- factor(0, levels=c(0,1))
rainy <- Mutasa5_new$rainy[!duplicated(Mutasa5_new$partid)]
parts <- levels(as.factor(Mutasa5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Mutasa5_clust_trips_shorter2_no_home3$rainy[which(Mutasa5_clust_trips_shorter2_no_home3$partid == part)] <- rainy[ii]
}

Mutasa5_clust_trips_shorter2_home4$rainy <- factor(0, levels=c(0,1))
rainy <- Mutasa5_new$rainy[!duplicated(Mutasa5_new$partid)]
parts <- levels(as.factor(Mutasa5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Mutasa5_clust_trips_shorter2_home4$rainy[which(Mutasa5_clust_trips_shorter2_home4$partid == part)] <- rainy[ii]
}
#
#### FOR NOW, REMOVE PARTS WITH RAINY AND DRY DATA ####
dry <- which(colSums(table(Mutasa5_new$dry,Mutasa5_new$partid)) == table(Mutasa5_new$dry,Mutasa5_new$partid)[2,]) ## 72 parts with data only in dry only
rainy <- which(colSums(table(Mutasa5_new$rainy,Mutasa5_new$partid)) == table(Mutasa5_new$rainy,Mutasa5_new$partid)[2,]) ## 74 parts with data in rainy only
both <- which(!(1:nlevels(as.factor(Mutasa5_new$partid)) %in% c(dry,rainy))) ## 36 parts with data in both seasons
rem <- levels(as.factor(Mutasa5_new$partid))[both]

Mutasa5_rain_b <- Mutasa5_new[-which(Mutasa5_new$partid %in% rem),]
Mutasa5_rain_clusts_shorter2 <- Mutasa5_clusts_shorter2[-which(Mutasa5_clusts_shorter2$partid %in% rem),]
Mutasa5_rain_clusts_shorter2_home <- Mutasa5_clusts_shorter2_home[-which(Mutasa5_clusts_shorter2_home$partid %in% rem),]
Mutasa5_rain_clusts_shorter2_no_home <- Mutasa5_clusts_shorter2_no_home[-which(Mutasa5_clusts_shorter2_no_home$partid %in% rem),]
Mutasa5_rain_clust_trips_shorter2_no_home <- Mutasa5_clust_trips_shorter2_no_home[-which(Mutasa5_clust_trips_shorter2_no_home$partid %in% rem),]
Mutasa5_rain_clust_trips_shorter2_home <- Mutasa5_clust_trips_shorter2_home[-which(Mutasa5_clust_trips_shorter2_home$partid %in% rem),]
Mutasa5_rain_clust_trips_shorter2_home3 <- Mutasa5_clust_trips_shorter2_home3[-which(Mutasa5_clust_trips_shorter2_home3$partid %in% rem),]
Mutasa5_rain_clust_trips_shorter2_no_home3 <- Mutasa5_clust_trips_shorter2_no_home3[-which(Mutasa5_clust_trips_shorter2_no_home3$partid %in% rem),]
Mutasa5_rain_clust_trips_shorter2_home4 <- Mutasa5_clust_trips_shorter2_home4[-which(Mutasa5_clust_trips_shorter2_home4$partid %in% rem),]
#
###### Split by Season ######
### number of locations (without home) ####
Mutasa5_rain_clusts_shorter2_no_home_dry <- Mutasa5_rain_clusts_shorter2_no_home[which(Mutasa5_rain_clusts_shorter2_no_home$rainy == 0),]
Mutasa5_rain_clusts_shorter2_no_home_rainy <- Mutasa5_rain_clusts_shorter2_no_home[which(Mutasa5_rain_clusts_shorter2_no_home$rainy == 1),]
loc_counts_dry <- favstats(Mutasa5_rain_clusts_shorter2_no_home_dry$clust~Mutasa5_rain_clusts_shorter2_no_home_dry$partid)$n
loc_counts_rainy <- favstats(Mutasa5_rain_clusts_shorter2_no_home_rainy$clust~Mutasa5_rain_clusts_shorter2_no_home_rainy$partid)$n
only_home_dry <- nlevels(as.factor(Mutasa5_rain_clusts_shorter2$partid[which(Mutasa5_rain_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Mutasa5_rain_clusts_shorter2_no_home_dry$partid))
only_home_rainy <- nlevels(as.factor(Mutasa5_rain_clusts_shorter2$partid[which(Mutasa5_rain_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Mutasa5_rain_clusts_shorter2_no_home_rainy$partid))
loc_counts_dry2 <- c(loc_counts_dry,rep(0, only_home_dry))
loc_counts_rainy2 <- c(loc_counts_rainy,rep(0, only_home_rainy))
favstats(loc_counts_dry2)
# min   Q1 median     Q3  max     mean       sd  n missing
#  0  4      8 11  25 8.041667 5.585286 72       0
favstats(loc_counts_rainy2)
# min Q1 median Q3 max     mean       sd  n missing
#   0   3      6 10.75  23 7.472973 5.527572 74       0

locs_count2 <- as.data.frame(cbind(c(loc_counts_dry2, loc_counts_rainy2),
                                   c(rep("dry",length(loc_counts_dry2)),
                                     rep("rainy",length(loc_counts_rainy2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-2,27)) +
  scale_x_continuous(breaks=c(seq(0,26,5)))
#

kruskal.test(dists_km2_season$values,dists_km2_season$type)
# p = 0.39
### distance of locations (without home) ####
### TOTAL ###
favstats(Mutasa5_rain_clusts_shorter2_no_home_dry$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#  2.745604e-05 0.4678129 1.32559 3.154556 481.4493 18.60404 60.2773 579       0
favstats(Mutasa5_rain_clusts_shorter2_no_home_rainy$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 0.006409266 0.5123223 1.283132 4.21008 481.4354 16.21704 55.10712 553       0

### PER PART ###
medians_km_dry <- (favstats(Mutasa5_rain_clusts_shorter2_no_home_dry$hhdist_km_haversine~Mutasa5_rain_clusts_shorter2_no_home_dry$partid)$median)
medians_km_rainy <- (favstats(Mutasa5_rain_clusts_shorter2_no_home_rainy$hhdist_km_haversine~Mutasa5_rain_clusts_shorter2_no_home_rainy$partid)$median)
favstats(medians_km_dry) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
# 0.08078032 0.6155971 1.139898 1.995036 299.3132 10.10343 40.03509 70       0
favstats(medians_km_rainy) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.1207797 0.6839959 1.177662 2.745979 217.5585 10.06372 36.09253 71       0

## REMOVE OUTLIER ## none
medians_km_dry_short <- medians_km_dry[-which(medians_km_dry > 100)]
medians_km_rainy_short <- medians_km_rainy[-which(medians_km_rainy > 100)]

favstats(medians_km_dry_short) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.08078032 0.6114316 1.007062 1.828719 50.98935 2.832056 7.581431 67       0
favstats(medians_km_rainy_short) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.1207797 0.6666691 1.073866 2.580591 57.42812 4.161585 9.257908 69       0

dists_km2_season <- as.data.frame(cbind(c(medians_km_dry_short, medians_km_rainy_short),
                                        c(rep("dry",length(medians_km_dry_short)),
                                          rep("rainy",length(medians_km_rainy_short)))))
colnames(dists_km2_season) <- c("values","type")
dists_km2_season$values <- as.numeric(as.character(dists_km2_season$values))
dists_km2_season$type <- factor(dists_km2_season$type)

ggplot(data=dists_km2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=10) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-2.5,65))

kruskal.test(dists_km2_season$values,dists_km2_season$type)
# p = 0.39
#
#########
### distance of locations (without home) ####
Mutasa5_clusts_shorter2_no_home_1_16_rainy <- Mutasa5_clusts_shorter2_no_home_1_16[which(Mutasa5_clusts_shorter2_no_home_1_16$rainy == 1),]
Mutasa5_clusts_shorter2_no_home_17_34_rainy <- Mutasa5_clusts_shorter2_no_home_17_34[which(Mutasa5_clusts_shorter2_no_home_17_34$rainy == 1),]
Mutasa5_clusts_shorter2_no_home_35_54_rainy <- Mutasa5_clusts_shorter2_no_home_35_54[which(Mutasa5_clusts_shorter2_no_home_35_54$rainy == 1),]
Mutasa5_clusts_shorter2_no_home_55_up_rainy <- Mutasa5_clusts_shorter2_no_home_55_up[which(Mutasa5_clusts_shorter2_no_home_55_up$rainy == 1),]
Mutasa5_clusts_shorter2_no_home_1_16_dry <- Mutasa5_clusts_shorter2_no_home_1_16[which(Mutasa5_clusts_shorter2_no_home_1_16$rainy == 0),]
Mutasa5_clusts_shorter2_no_home_17_34_dry <- Mutasa5_clusts_shorter2_no_home_17_34[which(Mutasa5_clusts_shorter2_no_home_17_34$rainy == 0),]
Mutasa5_clusts_shorter2_no_home_35_54_dry <- Mutasa5_clusts_shorter2_no_home_35_54[which(Mutasa5_clusts_shorter2_no_home_35_54$rainy == 0),]
Mutasa5_clusts_shorter2_no_home_55_up_dry <- Mutasa5_clusts_shorter2_no_home_55_up[which(Mutasa5_clusts_shorter2_no_home_55_up$rainy == 0),]



medians_km_1_16_rainy <- (favstats(Mutasa5_clusts_shorter2_no_home_1_16_rainy$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_1_16_rainy$partid)$median)
medians_km_17_34_rainy <- (favstats(Mutasa5_clusts_shorter2_no_home_17_34_rainy$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_17_34_rainy$partid)$median)
medians_km_35_54_rainy <- (favstats(Mutasa5_clusts_shorter2_no_home_35_54_rainy$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_35_54_rainy$partid)$median)
medians_km_55_up_rainy <- (favstats(Mutasa5_clusts_shorter2_no_home_55_up_rainy$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_55_up_rainy$partid)$median)
medians_km_1_16_dry <- (favstats(Mutasa5_clusts_shorter2_no_home_1_16_dry$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_1_16_dry$partid)$median)
medians_km_17_34_dry <- (favstats(Mutasa5_clusts_shorter2_no_home_17_34_dry$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_17_34_dry$partid)$median)
medians_km_35_54_dry <- (favstats(Mutasa5_clusts_shorter2_no_home_35_54_dry$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_35_54_dry$partid)$median)
medians_km_55_up_dry <- (favstats(Mutasa5_clusts_shorter2_no_home_55_up_dry$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_55_up_dry$partid)$median)

dists_km2_age <- as.data.frame(cbind(c(medians_km_1_16_rainy, medians_km_17_34_rainy, medians_km_35_54_rainy, medians_km_55_up_rainy,
                                       medians_km_1_16_dry, medians_km_17_34_dry, medians_km_35_54_dry, medians_km_55_up_dry),
                                     c(rep("1-16",length(medians_km_1_16_rainy)),
                                       rep("17-34",length(medians_km_17_34_rainy)),
                                       rep("35-54",length(medians_km_35_54_rainy)),
                                       rep("55+",length(medians_km_55_up_rainy)),
                                       rep("1-16",length(medians_km_1_16_dry)),
                                       rep("17-34",length(medians_km_17_34_dry)),
                                       rep("35-54",length(medians_km_35_54_dry)),
                                       rep("55+",length(medians_km_55_up_dry))),
                                     c( rep("rainy", (length(medians_km_1_16_rainy) + length(medians_km_17_34_rainy) + length(medians_km_35_54_rainy) + length(medians_km_55_up_rainy))),
                                        rep("dry", (length(medians_km_1_16_dry) + length(medians_km_17_34_dry) + length(medians_km_35_54_dry) + length(medians_km_55_up_dry))))))


colnames(dists_km2_age) <- c("values","type","rainy")
dists_km2_age$values <- as.numeric(as.character(dists_km2_age$values))
dists_km2_age$type <- factor(dists_km2_age$type)
dists_km2_age$rainy <- factor(dists_km2_age$rainy)

favstats(dists_km2_age$values[which(dists_km2_age$rainy=="rainy")]~dists_km2_age$type[which(dists_km2_age$rainy=="rainy")])
favstats(dists_km2_age$values[which(dists_km2_age$rainy=="dry")]~dists_km2_age$type[which(dists_km2_age$rainy=="dry")])

dists_km2_age_rainy <- dists_km2_age[which(dists_km2_age$rainy == "rainy"),]
dists_km2_age_dry <- dists_km2_age[which(dists_km2_age$rainy == "dry"),]


kruskal.test(dists_km2_age_rainy$values, dists_km2_age_rainy$type) # p = 0.23
kruskal.test(dists_km2_age_dry$values, dists_km2_age_dry$type) # p = 0.23


pairwise.wilcox.test(dists_km2_age_rainy$values, dists_km2_age_rainy$type, p.adjust.method = "holm")
pairwise.wilcox.test(dists_km2_age_dry$values, dists_km2_age_rainy$type, p.adjust.method = "holm")

cc <- rbind(dists_km2_age_dry, dists_km2_age_rainy)
kruskal.test(cc$values, cc$type) # p = 0.23
pairwise.wilcox.test(cc$values, cc$type, p.adjust.method = "holm")


ggplot(data=dists_km2_age_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=12) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-3,65))


ggplot(data=dists_km2_age_dry, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=12) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-3,65))
###
#########
### number of trips (without home) ####
### TOTAL ###
Mutasa5_rain_clust_trips_shorter2_no_home3_dry <- Mutasa5_rain_clust_trips_shorter2_no_home3[which(Mutasa5_rain_clust_trips_shorter2_no_home3$rainy == 0),]
Mutasa5_rain_clust_trips_shorter2_no_home3_rainy <- Mutasa5_rain_clust_trips_shorter2_no_home3[which(Mutasa5_rain_clust_trips_shorter2_no_home3$rainy == 1),]

favstats(Mutasa5_rain_clust_trips_shorter2_no_home3_dry$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      2  3  43 2.709845 3.929582 579       0
favstats(Mutasa5_rain_clust_trips_shorter2_no_home3_rainy$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      1  2  39 2.45208 3.313067 553       0

### PER PART ###
trips_part_median_dry <- favstats(Mutasa5_rain_clust_trips_shorter2_no_home3_dry$trip_count_new ~ Mutasa5_rain_clust_trips_shorter2_no_home3_dry$partid)$median
trips_part_median_rainy <- favstats(Mutasa5_rain_clust_trips_shorter2_no_home3_rainy$trip_count_new ~ Mutasa5_rain_clust_trips_shorter2_no_home3_rainy$partid)$median
favstats(trips_part_median_dry)
#    min Q1 median  Q3   max      mean       sd  n missing
#  1  1   1.75  2  13 1.885714 1.595218 70       0
favstats(trips_part_median_rainy)
#   min Q1 median  Q3   max      mean       sd  n missing
#  1  1    1.5  2  13 2.028169 2.028172 71       0

trips_part2_season <- as.data.frame(cbind(c(trips_part_median_dry, trips_part_median_rainy),
                                          c(rep("dry",length(trips_part_median_dry)),
                                            rep("rainy",length(trips_part_median_rainy)))))
colnames(trips_part2_season) <- c("values","type")
trips_part2_season$values <- as.numeric(as.character(trips_part2_season$values))
trips_part2_season$type <- factor(trips_part2_season$type)
kruskal.test(trips_part2_season$values,trips_part2_season$type)

ggplot(data=trips_part2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,15)) +
  scale_x_continuous(breaks=c(0:15))
#
### time per trip ####
### TOTAL ###
Mutasa5_rain_clust_trips_shorter2_no_home_dry <- Mutasa5_rain_clust_trips_shorter2_no_home[which(Mutasa5_rain_clust_trips_shorter2_no_home$rainy == 0),]
Mutasa5_rain_clust_trips_shorter2_no_home_rainy <- Mutasa5_rain_clust_trips_shorter2_no_home[which(Mutasa5_rain_clust_trips_shorter2_no_home$rainy == 1),]
favstats(Mutasa5_rain_clust_trips_shorter2_no_home_dry$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.48 1.436389 2.539444 5.445833 331.2 6.003608 17.73307 1569       0
favstats(Mutasa5_rain_clust_trips_shorter2_no_home_rainy$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
#  0.4413889 1.361597 2.661528 5.21 860.0692 7.311999 40.80837 1356       0

### PER PART ###
trip_time_part_median_dry <- favstats(Mutasa5_rain_clust_trips_shorter2_no_home_dry$trip_time_new_hour ~ Mutasa5_rain_clust_trips_shorter2_no_home_dry$partid)$median
trip_time_part_median_rainy <- favstats(Mutasa5_rain_clust_trips_shorter2_no_home_rainy$trip_time_new_hour ~ Mutasa5_rain_clust_trips_shorter2_no_home_rainy$partid)$median
favstats(trip_time_part_median_dry)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 1.146111 2.139896 2.783056 3.341215 33.96292 3.637617 4.024007 70       0
favstats(trip_time_part_median_rainy)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.066667 2.051806 3.066667 4.171389 12.13694 3.72552 2.394141 71       0

## REMOVE OUTLIERS ##
trip_time_part_median_dry_short <- trip_time_part_median_dry[which(trip_time_part_median_dry < 30)]
favstats(trip_time_part_median_dry_short)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.146111 2.133194 2.782778 3.287778 8.066389 3.19812 1.646474 69       0

trip_times_part2_season <- as.data.frame(cbind(c(trip_time_part_median_dry_short, trip_time_part_median_rainy),
                                               c(rep("dry",length(trip_time_part_median_dry_short)),
                                                 rep("rainy",length(trip_time_part_median_rainy)))))
colnames(trip_times_part2_season) <- c("values","type")
trip_times_part2_season$values <- as.numeric(as.character(trip_times_part2_season$values))
trip_times_part2_season$type <- factor(trip_times_part2_season$type)
≈
ggplot(data=trip_times_part2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,13)) +
  scale_x_continuous(breaks=c(seq(0,14,by=2)))
#
### Home Loc ####
Mutasa5_rain_clusts_shorter2_home_dry <- Mutasa5_rain_clusts_shorter2_home[which(Mutasa5_rain_clusts_shorter2_home$rainy == 0),] 
Mutasa5_rain_clusts_shorter2_home_rainy <- Mutasa5_rain_clusts_shorter2_home[which(Mutasa5_rain_clusts_shorter2_home$rainy == 1),] 

no_home_dry <- nlevels(as.factor(Mutasa5_rain_clusts_shorter2$partid[which(Mutasa5_rain_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Mutasa5_rain_clusts_shorter2_home_dry$partid))
## no drys with no home
no_home_rainy <- nlevels(as.factor(Mutasa5_rain_clusts_shorter2$partid[which(Mutasa5_rain_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Mutasa5_rain_clusts_shorter2_home_rainy$partid))
## no rainys with no home

### number of 'trips' to home loc ####
### TOTAL ###
Mutasa5_rain_clust_trips_shorter2_home3_dry <- Mutasa5_rain_clust_trips_shorter2_home3[which(Mutasa5_rain_clust_trips_shorter2_home3$rainy == 0),]
Mutasa5_rain_clust_trips_shorter2_home3_rainy <- Mutasa5_rain_clust_trips_shorter2_home3[which(Mutasa5_rain_clust_trips_shorter2_home3$rainy == 1),]

favstats(Mutasa5_rain_clust_trips_shorter2_home3_dry$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#   1  5     10 21  42 13.55172 10.3213 87       0
favstats(Mutasa5_rain_clust_trips_shorter2_home3_rainy$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#   1  5      8 17.75  35 11.60638 8.886643 94       0

### total time at home loc ####
Mutasa5_rain_clust_trips_shorter2_home4_dry <- Mutasa5_rain_clust_trips_shorter2_home4[which(Mutasa5_rain_clust_trips_shorter2_home4$rainy == 0),]
Mutasa5_rain_clust_trips_shorter2_home4_rainy <- Mutasa5_rain_clust_trips_shorter2_home4[which(Mutasa5_rain_clust_trips_shorter2_home4$rainy == 1),]

favstats(Mutasa5_rain_clust_trips_shorter2_home4_dry$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
# 0.920081 17.92198 23.60225 28.15999 65.96275 24.34449 11.49201 87       0
favstats(Mutasa5_rain_clust_trips_shorter2_home4_rainy$clust_trips_time_new)
#         min       Q1   median       Q3      max     mean       sd  n missing
# 6.550237 18.02888 25.23279 32.02426 60.75051 26.07171 11.07215 94       0
#
### percent time at home ####
favstats(Mutasa5_rain_clusts_shorter2_home_dry$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
#  29.7829 74.61344 87.68246 96.50758 100 81.55836 18.12059 72       0
## no one with <1% time at home (1 with <5%)
favstats(Mutasa5_rain_clusts_shorter2_home_rainy$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#  14.3927 73.73799 89.98518 96.17372 100 83.18846 18.26162 74       0
## no one with <5% time at home

perc_time_home2_season <- as.data.frame(cbind(c(Mutasa5_rain_clusts_shorter2_home_dry$perc_time_home, Mutasa5_rain_clusts_shorter2_home_rainy$perc_time_home),
                                              c(rep("dry",length(Mutasa5_rain_clusts_shorter2_home_dry$perc_time_home)),
                                                rep("rainy",length(Mutasa5_rain_clusts_shorter2_home_rainy$perc_time_home)))))
colnames(perc_time_home2_season) <- c("values","type")
perc_time_home2_season$values <- as.numeric(as.character(perc_time_home2_season$values))
perc_time_home2_season$type <- factor(perc_time_home2_season$type)

ggplot(data=perc_time_home2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,101))

kruskal.test(perc_time_home2_season$values,perc_time_home2_season$type)
# p = 0.59

##
### Average Direction of Travel #####
Mutasa5_rain_clusts_shorter2_no_home$rainy <- as.factor(Mutasa5_rain_clusts_shorter2_no_home$rainy)
temp <- Mutasa5_rain_clusts_shorter2_no_home[!duplicated(Mutasa5_rain_clusts_shorter2_no_home$partid),]
temp2 <- temp[which(!(is.na(temp$rainy))),]

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=rainy), position=position_dodge2(preserve = "single"), width=4) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Rainy Season", values=c("sienna","skyblue"))
#

favstats(temp2$avg_dir_from_home_bearing~temp2$rainy)
#   temp2$rainy      min       Q1   median       Q3      max     mean        sd  n missing
# 1           0 2.636555 127.3184 229.1859 288.7319 359.1797 208.9773 100.81737 86       0
# 2           1 6.865241 112.5242 217.5350 264.0063 359.7089 196.8002  99.29785 90       0
favstats(temp2$avg_dir_from_home_mag~temp2$rainy)
#   temp2$rainy        min        Q1   median       Q3      max     mean       sd  n missing
# 1           0 0.03380367 0.5417520 1.538976 6.991581 210.0080 16.84851 39.23349 86       0
# 2           1 0.08965062 0.7174819 1.303399 3.756298 186.6763 10.10687 29.50991 90       0

temp3 <- temp2[which(temp2$avg_dir_from_home_mag < 70),] ## 164 of 176 parts

favstats(temp3$avg_dir_from_home_bearing~temp3$rainy)
#   temp3$rainy      min       Q1   median       Q3      max     mean       sd  n missing
# 1           0 2.636555 119.4387 221.1990 292.3151 359.1797 204.3677 101.7149 77       0
# 2           1 6.865241 106.8027 211.7342 262.5723 359.7089 194.4994 100.0912 87       0
favstats(temp3$avg_dir_from_home_mag~temp3$rainy)
#   temp3$rainy        min        Q1   median       Q3      max     mean       sd  n missing
# 1           0 0.03380367 0.5265798 1.244584 4.026456 65.41056 4.983014 11.09078 77       0
# 2           1 0.08965062 0.6733483 1.276329 3.537346 68.50646 5.305079 10.77333 87       0

kruskal.test(temp3$avg_dir_from_home_bearing,temp3$rainy)

ggplot(data=temp3, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=rainy), position=position_dodge2(preserve = "single"), width=4) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Rainy Season", values=c("sienna","skyblue"))




temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)] <- -(360-temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)])
ggplot(data=temp2, aes(x=avg_dir_from_home_values, y=..count.., fill=rainy)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~rainy, labeller = labeller(rainy = c("0"="Dry", "1"="Rainy")))+
  scale_fill_manual(values=c("sienna","skyblue")) 

#

temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values < 0)] <- (360+temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values < 0)])

## average direction calc
temp2$values2 <- 90 - temp2$avg_dir_from_home_values
temp2$values2[which(temp2$values2 < 0)] <- temp2$values2[which(temp2$values2 < 0)] + 360
temp2$values3 <- deg2rad(temp2$values2)

temp2$avg_value <- as.numeric(NA)
for(ii in 1:nlevels(temp2$rainy)){
  level <- levels(temp2$rainy)[ii]
  temp <- temp2[which(temp2$rainy == level),]
  sum_cos_dist <- sum(cos(temp$values3))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$values3))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  temp2$avg_value[which(temp2$rainy == level)] <- avg_deg
}

temp2$avg_value2 <- 90 - temp2$avg_value
temp2$avg_value2[which(temp2$avg_value2 < 0)] <- temp2$avg_value2[which(temp2$avg_value2 < 0)] + 360
favstats(temp2$avg_value2~temp2$rainy)
## dry: 305.0 --> NW
## rainy: 267.4 --> W



#
### time spent at clusters vs elsewhere (travel, less important places, etc) ####
# percent time in clusters
Mutasa5_rain_clusts_shorter2_dry <- Mutasa5_rain_clusts_shorter2[which(Mutasa5_rain_clusts_shorter2$rainy == 0),]
Mutasa5_rain_clusts_shorter2_rainy <- Mutasa5_rain_clusts_shorter2[which(Mutasa5_rain_clusts_shorter2$rainy == 1),]
Mutasa5_rain_clusts_shorter2_dry$perc_time_in_all_clusts <- Mutasa5_rain_clusts_shorter2_dry$part_trips_time_new/Mutasa5_rain_clusts_shorter2_dry$total_time
Mutasa5_rain_clusts_shorter2_rainy$perc_time_in_all_clusts <- Mutasa5_rain_clusts_shorter2_rainy$part_trips_time_new/Mutasa5_rain_clusts_shorter2_rainy$total_time

summary(Mutasa5_rain_clusts_shorter2_dry$perc_time_in_all_clusts[!duplicated(Mutasa5_rain_clusts_shorter2_dry$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6629  0.8500  0.9261  0.8949  0.9534  0.9955 
summary(Mutasa5_rain_clusts_shorter2_rainy$perc_time_in_all_clusts[!duplicated(Mutasa5_rain_clusts_shorter2_rainy$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8287  0.8902  0.9373  0.9318  0.9765  0.9950 

# time outside clusters
summary(Mutasa5_rain_clusts_shorter2_dry$total_time[!duplicated(Mutasa5_rain_clusts_shorter2_dry$partid)] - Mutasa5_rain_clusts_shorter2_dry$part_trips_time_new[!duplicated(Mutasa5_rain_clusts_shorter2_dry$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1395  1.1869  2.2344  3.2449  4.6457  9.8070 
summary(Mutasa5_rain_clusts_shorter2_rainy$total_time[!duplicated(Mutasa5_rain_clusts_shorter2_rainy$partid)] - Mutasa5_rain_clusts_shorter2_rainy$part_trips_time_new[!duplicated(Mutasa5_rain_clusts_shorter2_rainy$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2383  0.8182  1.8176  2.4824  3.9561  6.6361 
#
#
#
#
### Time at locations ####
favstats(Mutasa5_rain_clusts_shorter2_no_home_dry$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.5072222 2.234028 4.173889 9.701944 570.4803 16.26884 47.81858 579       0
favstats(Mutasa5_rain_clusts_shorter2_no_home_rainy$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.5241667 2.310833 4.021667 8.451667 1537.291 17.9296 86.66931 553       0

### PER PART ###
clust_time_part_median_dry <- (favstats(Mutasa5_rain_clusts_shorter2_no_home_dry$clust_trips_time_new~Mutasa5_rain_clusts_shorter2_no_home_dry$partid)$median)*24
clust_time_part_median_rainy <- (favstats(Mutasa5_rain_clusts_shorter2_no_home_rainy$clust_trips_time_new~Mutasa5_rain_clusts_shorter2_no_home_rainy$partid)$median)*24
favstats(clust_time_part_median_dry)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 2.020556 3.294097 4.376944 6.006354 195.3874 9.89598 26.78139 70       0
favstats(clust_time_part_median_rainy)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.656389 3.257778 4.142222 6.609167 98.0725 8.882533 16.26301 71       0
clust_time_part_median_dry2 <- clust_time_part_median_dry[which(clust_time_part_median_dry < 150)]
favstats(clust_time_part_median_dry2)
# min       Q1   median Q3      max     mean       sd  n missing
# 2.020556 3.283889 4.354861  6 110.4386 7.207699 14.64485 69       0

clust_time_part_rainy <- as.data.frame(cbind(c(clust_time_part_median_dry2, clust_time_part_median_rainy),
                                             c(rep("dry",length(clust_time_part_median_dry2)),
                                               rep("rainy",length(clust_time_part_median_rainy)))))
colnames(clust_time_part_rainy) <- c("values","type")
clust_time_part_rainy$values <- as.numeric(as.character(clust_time_part_rainy$values))
clust_time_part_rainy$type <- factor(clust_time_part_rainy$type)
kruskal.test(clust_time_part_rainy$values,clust_time_part_rainy$type)

ggplot(data=clust_time_part_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=5) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-1,110)) +
  scale_x_continuous(breaks=c(seq(0,120,10)))
#
#
#
####### Biting time metrics #####
### number locations outside home ####
Mutasa5_clusts_shorter2_no_home_biting_places <- Mutasa5_clusts_shorter2_no_home[which(!(is.na(Mutasa5_clusts_shorter2_no_home$clust_trips_time_new_biting_time))),]
Mutasa5_clusts_shorter2_no_home_biting_places2 <- Mutasa5_clusts_shorter2_no_home_biting_places[which((Mutasa5_clusts_shorter2_no_home_biting_places$clust_trips_time_new_biting_time*24) >= 0.5),]

Mutasa5_clusts_shorter2_no_home_biting_places2_dry <- Mutasa5_clusts_shorter2_no_home_biting_places2[which(Mutasa5_clusts_shorter2_no_home_biting_places2$rainy == 0),]
Mutasa5_clusts_shorter2_no_home_biting_places2_rainy <- Mutasa5_clusts_shorter2_no_home_biting_places2[which(Mutasa5_clusts_shorter2_no_home_biting_places2$rainy == 1),]

loc_counts_dry <- favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_dry$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_dry$partid)$n
loc_counts_rainy <- favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_rainy$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_rainy$partid)$n
only_home_dry <- nlevels(as.factor(Mutasa5_rain_clusts_shorter2$partid[which(Mutasa5_rain_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_dry$partid))
only_home_rainy <- nlevels(as.factor(Mutasa5_rain_clusts_shorter2$partid[which(Mutasa5_rain_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_rainy$partid))
loc_counts_dry2 <- c(loc_counts_dry,rep(0, only_home_dry))
loc_counts_rainy2 <- c(loc_counts_rainy,rep(0, only_home_rainy))

favstats(loc_counts_dry2)
# min Q1 median Q3 max     mean       sd  n missing
#  0 0.75      1  2   7 1.681 1.617 72       0
favstats(loc_counts_rainy2)
# min Q1 median Q3 max     mean       sd  n missing
#    0  0      1  2   8 1.784 1.995 74       0
summary(as.factor(loc_counts_dry2))
#  0  1  2  3  4  5  6  7 
# 18 22 16  5  7  1  2  1 
summary(as.factor(loc_counts_rainy2))
#  0  1  2  3  4  6  8 
# 20 21 16  9  2  2  4 
## all those with 0 have biting time at home location (not just nowhere)

locs_count2 <- as.data.frame(cbind(c(loc_counts_dry2, loc_counts_rainy2),
                                   c(rep("dry",length(loc_counts_dry2)),
                                     rep("rainy",length(loc_counts_rainy2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-0.35,8.35)) +
  scale_x_continuous(breaks=c(seq(0,8,1)))
#
##
### time per location ####
favstats((Mutasa5_clusts_shorter2_no_home_biting_places2_dry$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5274 1.638   7.92 15.57 196.4 18.68 33.37 121       0
favstats((Mutasa5_clusts_shorter2_no_home_biting_places2_rainy$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5072 2.116  7.228 15.56 504.1 23.61 58.08 132       0

clusts_median_Mutasa_dry <- (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_dry$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_dry$partid)$median)*24
clusts_median_Mutasa_rainy <- (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_rainy$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_rainy$partid)$median)*24
favstats(clusts_median_Mutasa_dry)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.5831 2.493  7.991 12.99 130.5   15 24.4 54       0
favstats(clusts_median_Mutasa_rainy)
# min       Q1   median  Q3      max     mean       sd  n missing
# 1.012 3.88  8.614 32.35 504.1 34.06 78.18 54       0

clust_time_part_rainy <- as.data.frame(cbind(c(clusts_median_Mutasa_dry, clusts_median_Mutasa_rainy),
                                             c(rep("dry",length(clusts_median_Mutasa_dry)),
                                               rep("rainy",length(clusts_median_Mutasa_rainy)))))
colnames(clust_time_part_rainy) <- c("values","type")
clust_time_part_rainy$values <- as.numeric(as.character(clust_time_part_rainy$values))
clust_time_part_rainy$type <- factor(clust_time_part_rainy$type)

ggplot(data=clust_time_part_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=2.5) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,500)) +
  scale_x_continuous(breaks=c(seq(0,500,50)))

clust_time_part_rainy2 <- clust_time_part_rainy[which(clust_time_part_rainy$values < 500),]
ggplot(data=clust_time_part_rainy2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=2.5) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,250)) +
  scale_x_continuous(breaks=c(seq(0,250,50)))
#
### number trips ####
Mutasa5_clust_trips_shorter2_no_home_biting_places <- Mutasa5_clust_trips_shorter2_no_home[which(!(is.na(Mutasa5_clust_trips_shorter2_no_home$trip_time_new2_biting_time))),]
## only when >30 minutes during biting time
Mutasa5_clust_trips_shorter2_no_home_biting_places2 <- Mutasa5_clust_trips_shorter2_no_home_biting_places[which((Mutasa5_clust_trips_shorter2_no_home_biting_places$trip_time_new2_biting_time*24) >= 0.5),]
## removes 63 trips, 3 participants


Mutasa5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2 <- numeric(length=nrow(Mutasa5_clust_trips_shorter2_no_home_biting_places2))
parts <- levels(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places2$clust[which(Mutasa5_clust_trips_shorter2_no_home_biting_places2$partid == part)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    counts <- Mutasa5_clust_trips_shorter2_no_home_biting_places2$trip_count_new[which(Mutasa5_clust_trips_shorter2_no_home_biting_places2$partid == part & Mutasa5_clust_trips_shorter2_no_home_biting_places2$clust == clust)]
    counts <- as.factor(counts)
    levels(counts) <- c(1:nlevels(counts))
    Mutasa5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2[which(Mutasa5_clust_trips_shorter2_no_home_biting_places2$partid == part & Mutasa5_clust_trips_shorter2_no_home_biting_places2$clust == clust)] <- as.numeric(as.character(counts))
  }
}

Mutasa5_clust_trips_shorter2_no_home_biting_places2b <- Mutasa5_clust_trips_shorter2_no_home_biting_places2[order(Mutasa5_clust_trips_shorter2_no_home_biting_places2$partid,
                                                                                                                Mutasa5_clust_trips_shorter2_no_home_biting_places2$clust,
                                                                                                                -Mutasa5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2),]
Mutasa5_clust_trips_shorter2_no_home_biting_places3 <- Mutasa5_clust_trips_shorter2_no_home_biting_places2b[!duplicated(Mutasa5_clust_trips_shorter2_no_home_biting_places2b[,c(1,15)]),]

Mutasa5_clust_trips_shorter2_no_home_biting_places3_dry <- Mutasa5_clust_trips_shorter2_no_home_biting_places3[which(Mutasa5_clust_trips_shorter2_no_home_biting_places3$rainy == 0),]
Mutasa5_clust_trips_shorter2_no_home_biting_places3_rainy <- Mutasa5_clust_trips_shorter2_no_home_biting_places3[which(Mutasa5_clust_trips_shorter2_no_home_biting_places3$rainy == 1),]

favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2  22 2.124 3.081 121       0
favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#  1  1      1  2  22 2.25 3.204 132       0

summary(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2))
#   1  2  3  4  6  7  8  9 10 13 14 22 
# 87 17  3  4  1  1  2  1  1  1  2  1 
summary(as.factor(Mutasa5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2))
#  1  2  3  4  5  6  7  8  9 11 22 
# 88 19 10  2  1  2  1  2  3  2  2 

### PER PART ###
trip_count_median_dry <- favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_dry$partid)$median
favstats(trip_count_median_dry)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1 1.5  13 1.806 2.245 54       0
trip_count_median_rainy <- favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_rainy$partid)$median
favstats(trip_count_median_rainy)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2  11 1.796 1.67 54       0

trips_part2_rainy <- as.data.frame(cbind(c(trip_count_median_dry, trip_count_median_rainy),
                                         c(rep("dry",length(trip_count_median_dry)),
                                           rep("rainy",length(trip_count_median_rainy)))))
colnames(trips_part2_rainy) <- c("values","type")
trips_part2_rainy$values <- as.numeric(as.character(trips_part2_rainy$values))
trips_part2_rainy$type <- factor(trips_part2_rainy$type)

ggplot(data=trips_part2_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2.5) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,13.25)) + 
  scale_x_continuous(breaks=c(seq(0,14,2)))


## just medians
ggplot(data=trips_part2[which(trips_part2$type=="part_median"),], aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,30))+ 
  scale_x_continuous(breaks=c(seq(0,30,2)))
#
##### percent biting time spent at home vs elsewhere #####
Mutasa5_clusts_shorter2_biting_places <- Mutasa5_clusts_shorter2[which(!(is.na(Mutasa5_clusts_shorter2$clust_trips_time_new_biting_time))),]
Mutasa5_clusts_shorter2_biting_places$percent_home_biting_time <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_biting_places$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  percent <- Mutasa5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Mutasa5_clusts_shorter2_biting_places$partid == part & Mutasa5_clusts_shorter2_biting_places$home_clust2 == 1)]/
    sum(Mutasa5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Mutasa5_clusts_shorter2_biting_places$partid == part)])
  if(length(percent) == 0){
    percent <- 0
  }
  Mutasa5_clusts_shorter2_biting_places$percent_home_biting_time[which(Mutasa5_clusts_shorter2_biting_places$partid == part)] <- percent*100
}
Mutasa5_clusts_shorter2_biting_places2 <- Mutasa5_clusts_shorter2_biting_places[!duplicated(Mutasa5_clusts_shorter2_biting_places$partid),]
## none with 0%
homes <- Mutasa5_clusts_shorter2_home$partid
Mutasa5_clusts_shorter2_biting_places3 <- Mutasa5_clusts_shorter2_biting_places[which(Mutasa5_clusts_shorter2_biting_places$partid %in% homes),]
Mutasa5_clusts_shorter2_biting_places3b <- Mutasa5_clusts_shorter2_biting_places3[!duplicated(Mutasa5_clusts_shorter2_biting_places3$partid),]


Mutasa5_clusts_shorter2_biting_places3b_dry <- Mutasa5_clusts_shorter2_biting_places3b[which(Mutasa5_clusts_shorter2_biting_places3b$rainy == 0),]
Mutasa5_clusts_shorter2_biting_places3b_rainy <- Mutasa5_clusts_shorter2_biting_places3b[which(Mutasa5_clusts_shorter2_biting_places3b$rainy == 1),]

favstats(Mutasa5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time)
favstats(Mutasa5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time[which(Mutasa5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time != 0)])
favstats(Mutasa5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time)
favstats(Mutasa5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time[which(Mutasa5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time != 0)])


perc_time_home2_rainy <- as.data.frame(cbind(c(Mutasa5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time, Mutasa5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time),
                                             c(rep("dry",length(Mutasa5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time)),
                                               rep("rainy",length(Mutasa5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time)))))
colnames(perc_time_home2_rainy) <- c("values","type")
perc_time_home2_rainy$values <- as.numeric(as.character(perc_time_home2_rainy$values))
perc_time_home2_rainy$type <- factor(perc_time_home2_rainy$type)

ggplot(data=perc_time_home2_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-1,101))
##
###############
###############
##### NCHELENGE #####
Nchelenge5_new$rainy <- factor(0, levels=c(0,1))
Nchelenge5_new$dry <- factor(0, levels=c(0,1))
Nchelenge5_new$month <- as.numeric(c(unlist(strsplit(Nchelenge5_new$date_new,"-"))[c(seq(from=2,to=c((nrow(Nchelenge5_new)*3)-1), by=3))]))
Nchelenge5_new$rainy[which(Nchelenge5_new$month %in% c(1,2,3,10,11,12))] <- 1
Nchelenge5_new$dry[which(Nchelenge5_new$month %in% c(4:9))] <- 1
table(Nchelenge5_new$dry,Nchelenge5_new$partid) ## 15287 is only part with data in both rainy and dry
dry <- which(colSums(table(Nchelenge5_new$dry,Nchelenge5_new$partid)) == table(Nchelenge5_new$dry,Nchelenge5_new$partid)[2,]) ## 37 parts with data only in 1 season
rainy <- which(colSums(table(Nchelenge5_new$rainy,Nchelenge5_new$partid)) == table(Nchelenge5_new$rainy,Nchelenge5_new$partid)[2,]) ## 38 parts with data in both seasons
both <- which(!(1:nlevels(as.factor(Nchelenge5_new$partid)) %in% c(dry,rainy))) ## no parts with data in both seasons

Nchelenge5_clusts_shorter2$rainy <- factor(0, levels=c(0,1))
rainy <- Nchelenge5_new$rainy[!duplicated(Nchelenge5_new$partid)]
parts <- levels(as.factor(Nchelenge5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Nchelenge5_clusts_shorter2$rainy[which(Nchelenge5_clusts_shorter2$partid == part)] <- rainy[ii]
}

Nchelenge5_clusts_shorter2_home$rainy <- factor(0, levels=c(0,1))
rainy <- Nchelenge5_new$rainy[!duplicated(Nchelenge5_new$partid)]
parts <- levels(as.factor(Nchelenge5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Nchelenge5_clusts_shorter2_home$rainy[which(Nchelenge5_clusts_shorter2_home$partid == part)] <- rainy[ii]
}

Nchelenge5_clusts_shorter2_no_home$rainy <- factor(0, levels=c(0,1))
rainy <- Nchelenge5_new$rainy[!duplicated(Nchelenge5_new$partid)]
parts <- levels(as.factor(Nchelenge5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Nchelenge5_clusts_shorter2_no_home$rainy[which(Nchelenge5_clusts_shorter2_no_home$partid == part)] <- rainy[ii]
}

Nchelenge5_clust_trips_shorter2_no_home$rainy <- factor(0, levels=c(0,1))
rainy <- Nchelenge5_new$rainy[!duplicated(Nchelenge5_new$partid)]
parts <- levels(as.factor(Nchelenge5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Nchelenge5_clust_trips_shorter2_no_home$rainy[which(Nchelenge5_clust_trips_shorter2_no_home$partid == part)] <- rainy[ii]
}

Nchelenge5_clust_trips_shorter2_home$rainy <- factor(0, levels=c(0,1))
rainy <- Nchelenge5_new$rainy[!duplicated(Nchelenge5_new$partid)]
parts <- levels(as.factor(Nchelenge5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Nchelenge5_clust_trips_shorter2_home$rainy[which(Nchelenge5_clust_trips_shorter2_home$partid == part)] <- rainy[ii]
}

Nchelenge5_clust_trips_shorter2_home3$rainy <- factor(0, levels=c(0,1))
rainy <- Nchelenge5_new$rainy[!duplicated(Nchelenge5_new$partid)]
parts <- levels(as.factor(Nchelenge5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Nchelenge5_clust_trips_shorter2_home3$rainy[which(Nchelenge5_clust_trips_shorter2_home3$partid == part)] <- rainy[ii]
}

Nchelenge5_clust_trips_shorter2_no_home3$rainy <- factor(0, levels=c(0,1))
rainy <- Nchelenge5_new$rainy[!duplicated(Nchelenge5_new$partid)]
parts <- levels(as.factor(Nchelenge5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Nchelenge5_clust_trips_shorter2_no_home3$rainy[which(Nchelenge5_clust_trips_shorter2_no_home3$partid == part)] <- rainy[ii]
}

Nchelenge5_clust_trips_shorter2_home4$rainy <- factor(0, levels=c(0,1))
rainy <- Nchelenge5_new$rainy[!duplicated(Nchelenge5_new$partid)]
parts <- levels(as.factor(Nchelenge5_new$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  Nchelenge5_clust_trips_shorter2_home4$rainy[which(Nchelenge5_clust_trips_shorter2_home4$partid == part)] <- rainy[ii]
}
#
###### Split by Season ######
### number of locations (without home) ####
Nchelenge5_clusts_shorter2_no_home_dry <- Nchelenge5_clusts_shorter2_no_home[which(Nchelenge5_clusts_shorter2_no_home$rainy == 0),]
Nchelenge5_clusts_shorter2_no_home_rainy <- Nchelenge5_clusts_shorter2_no_home[which(Nchelenge5_clusts_shorter2_no_home$rainy == 1),]
loc_counts_dry <- favstats(Nchelenge5_clusts_shorter2_no_home_dry$clust~Nchelenge5_clusts_shorter2_no_home_dry$partid)$n
loc_counts_rainy <- favstats(Nchelenge5_clusts_shorter2_no_home_rainy$clust~Nchelenge5_clusts_shorter2_no_home_rainy$partid)$n
only_home_dry <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_dry$partid))
only_home_rainy <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_rainy$partid))
loc_counts_dry2 <- c(loc_counts_dry,rep(0, only_home_dry))
loc_counts_rainy2 <- c(loc_counts_rainy,rep(0, only_home_rainy))
favstats(loc_counts_dry2)
# min   Q1 median     Q3  max     mean       sd  n missing
#   0  5      9 15  23 9.648649 5.99637 37       0
favstats(loc_counts_rainy2)
# min Q1 median Q3 max     mean       sd  n missing
#   0  3    7.5 11  19 7.526316 5.23882 38       0

locs_count2 <- as.data.frame(cbind(c(loc_counts_dry2, loc_counts_rainy2),
                                   c(rep("dry",length(loc_counts_dry2)),
                                     rep("rainy",length(loc_counts_rainy2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-2,26)) +
  scale_x_continuous(breaks=c(seq(0,26,5)))
#
kruskal.test(locs_count2$values,locs_count2$type)
# p=0.133

### distance of locations (without home) ####
### TOTAL ###
favstats(Nchelenge5_clusts_shorter2_no_home_dry$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
# 0.02460256 0.3409836 1.383268 4.792105 187.2311 8.385091 25.51193 357       0
favstats(Nchelenge5_clusts_shorter2_no_home_rainy$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
#  0 0.3879079 1.169316 4.180639 55.26877 4.93061 9.683853 286       0

### PER PART ###
medians_km_dry <- (favstats(Nchelenge5_clusts_shorter2_no_home_dry$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_dry$partid)$median)
medians_km_rainy <- (favstats(Nchelenge5_clusts_shorter2_no_home_rainy$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_rainy$partid)$median)
favstats(medians_km_dry) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
# 0.1012507 0.539771 0.9653387 2.41848 186.7748 9.463793 33.66226 36       0
favstats(medians_km_rainy) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.1811803 0.5871372 1.385886 3.867751 25.71711 3.72294 5.667487 35       0

## REMOVE OUTLIER ##
medians_km_dry_short <- medians_km_dry[-which(medians_km_dry > 70)]
favstats(medians_km_dry_short)
# min        Q1    median       Q3      max     mean       sd  n missing
# 0.1012507 0.5376021 0.9457531 2.053471 24.49658 2.014665 4.136003 34       0

medians_km_rainy_short <- medians_km_rainy

dists_km2_season <- as.data.frame(cbind(c(medians_km_dry_short, medians_km_rainy_short),
                                        c(rep("dry",length(medians_km_dry_short)),
                                          rep("rainy",length(medians_km_rainy_short)))))
colnames(dists_km2_season) <- c("values","type")
dists_km2_season$values <- as.numeric(as.character(dists_km2_season$values))
dists_km2_season$type <- factor(dists_km2_season$type)

ggplot(data=dists_km2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,30))

kruskal.test(dists_km2_season$values,dists_km2_season$type)
# p=0.12




dists_km2_season2 <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_no_home_dry$hhdist_km_haversine, Nchelenge5_clusts_shorter2_no_home_rainy$hhdist_km_haversine),
                                        c(rep("dry",length(Nchelenge5_clusts_shorter2_no_home_dry$hhdist_km_haversine)),
                                          rep("rainy",length(Nchelenge5_clusts_shorter2_no_home_rainy$hhdist_km_haversine)))))
colnames(dists_km2_season2) <- c("values","type")
dists_km2_season2$values <- as.numeric(as.character(dists_km2_season2$values))
dists_km2_season2$type <- factor(dists_km2_season2$type)

ggplot(data=dists_km2_season2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,70))

#
### number of trips (without home) ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2_no_home3_dry <- Nchelenge5_clust_trips_shorter2_no_home3[which(Nchelenge5_clust_trips_shorter2_no_home3$rainy == 0),]
Nchelenge5_clust_trips_shorter2_no_home3_rainy <- Nchelenge5_clust_trips_shorter2_no_home3[which(Nchelenge5_clust_trips_shorter2_no_home3$rainy == 1),]

favstats(Nchelenge5_clust_trips_shorter2_no_home3_dry$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      2  3  23 2.662011 2.944921 358       0
favstats(Nchelenge5_clust_trips_shorter2_no_home3_rainy$trip_count_new)
# min Q1 median  Q3 max     mean       sd     n missing
#    1  1      1  3  17 2.332168 2.382277 286       0

### PER PART ###
trips_part_median_dry <- favstats(Nchelenge5_clust_trips_shorter2_no_home3_dry$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_dry$partid)$median
trips_part_median_rainy <- favstats(Nchelenge5_clust_trips_shorter2_no_home3_rainy$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_rainy$partid)$median
favstats(trips_part_median_dry)
#    min Q1 median  Q3   max      mean       sd  n missing
#  1  1    1.5  2   4 1.638889 0.6929349 36       0
favstats(trips_part_median_rainy)
#   min Q1 median  Q3   max      mean       sd  n missing
#    1  1      2  2   3 1.642857 0.536484 35       0

trips_part2_season <- as.data.frame(cbind(c(trips_part_median_dry, trips_part_median_rainy),
                                          c(rep("dry",length(trips_part_median_dry)),
                                            rep("rainy",length(trips_part_median_rainy)))))
colnames(trips_part2_season) <- c("values","type")
trips_part2_season$values <- as.numeric(as.character(trips_part2_season$values))
trips_part2_season$type <- factor(trips_part2_season$type)

ggplot(data=trips_part2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks=c(0:10))
#
kruskal.test(trips_part2_season$values,trips_part2_season$type)
# p=0.72
### time per trip ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2_no_home_dry <- Nchelenge5_clust_trips_shorter2_no_home[which(Nchelenge5_clust_trips_shorter2_no_home$rainy == 0),]
Nchelenge5_clust_trips_shorter2_no_home_rainy <- Nchelenge5_clust_trips_shorter2_no_home[which(Nchelenge5_clust_trips_shorter2_no_home$rainy == 1),]
favstats(Nchelenge5_clust_trips_shorter2_no_home_dry$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2416667 0.5519444 1.163056 2.584167 329.3114 5.10986 20.35984 945       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_rainy$trip_time_new_hour)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2466667 0.5941667 1.241944 2.902778 394.7653 6.4893 28.5097 665       0

### PER PART ###
trip_time_part_median_dry <- favstats(Nchelenge5_clust_trips_shorter2_no_home_dry$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_dry$partid)$median
trip_time_part_median_rainy <- favstats(Nchelenge5_clust_trips_shorter2_no_home_rainy$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_rainy$partid)$median
favstats(trip_time_part_median_dry)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4955556 0.8979861 1.188819 1.743229 6.148333 1.537959 1.062494 36       0
favstats(trip_time_part_median_rainy)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.4844444 0.878125 1.391806 1.869722 6.663194 1.621921 1.16963 35       0

## REMOVE OUTLIERS ##
trip_time_part_median_dry_short <- trip_time_part_median_dry ## no outliers to remove

trip_times_part2_season <- as.data.frame(cbind(c(trip_time_part_median_dry_short, trip_time_part_median_rainy),
                                               c(rep("dry",length(trip_time_part_median_dry_short)),
                                                 rep("rainy",length(trip_time_part_median_rainy)))))
colnames(trip_times_part2_season) <- c("values","type")
trip_times_part2_season$values <- as.numeric(as.character(trip_times_part2_season$values))
trip_times_part2_season$type <- factor(trip_times_part2_season$type)

ggplot(data=trip_times_part2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 0.5, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,7)) +
  scale_x_continuous(breaks=c(0:7))
#
kruskal.test(trip_times_part2_season$values,trip_times_part2_season$type)
# p =0.73
### Home Loc ####
Nchelenge5_clusts_shorter2_home_dry <- Nchelenge5_clusts_shorter2_home[which(Nchelenge5_clusts_shorter2_home$rainy == 0),] 
Nchelenge5_clusts_shorter2_home_rainy <- Nchelenge5_clusts_shorter2_home[which(Nchelenge5_clusts_shorter2_home$rainy == 1),] 

no_home_dry <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_home_dry$partid))
## no drys with no home
no_home_rainy <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_home_rainy$partid))
## no rainys with no home

### number of 'trips' to home loc ####
### TOTAL ###
Nchelenge5_clust_trips_shorter2_home3_dry <- Nchelenge5_clust_trips_shorter2_home3[which(Nchelenge5_clust_trips_shorter2_home3$rainy == 0),]
Nchelenge5_clust_trips_shorter2_home3_rainy <- Nchelenge5_clust_trips_shorter2_home3[which(Nchelenge5_clust_trips_shorter2_home3$rainy == 1),]

favstats(Nchelenge5_clust_trips_shorter2_home3_dry$trip_count_new)
#   min Q1 median Q3 max     mean       sd  n missing
#   1  5     10 15  48 12.02703 9.967967 37       0
favstats(Nchelenge5_clust_trips_shorter2_home3_rainy$trip_count_new)
#   min Q1 median    Q3  max     mean       sd  n missing
#   1  3      5  9  33 9.026316 9.086631 38       0

### total time at home loc ####
Nchelenge5_clust_trips_shorter2_home4_dry <- Nchelenge5_clust_trips_shorter2_home4[which(Nchelenge5_clust_trips_shorter2_home4$rainy == 0),]
Nchelenge5_clust_trips_shorter2_home4_rainy <- Nchelenge5_clust_trips_shorter2_home4[which(Nchelenge5_clust_trips_shorter2_home4$rainy == 1),]

favstats(Nchelenge5_clust_trips_shorter2_home4_dry$clust_trips_time_new)
#        min       Q1   median       Q3      max     mean       sd  n missing
# 0.414523 14.41278 22.18398 28.5879 33.95134 20.82611 9.294915 37       0
favstats(Nchelenge5_clust_trips_shorter2_home4_rainy$clust_trips_time_new)
#         min       Q1   median       Q3      max     mean       sd  n missing
# 0.028894 16.51015 24.39648 30.56072 40.73759 22.60721 10.71608 38       0
#
### percent time at home ####

favstats(Nchelenge5_clusts_shorter2_home_dry$perc_time_home)
#  min       Q1   median       Q3    max    mean       sd  n missing
#  1.618541 71.84745 89.15795 97.89361 100 80.36563 23.97399 37       0
favstats(Nchelenge5_clusts_shorter2_home_rainy$perc_time_home)
#  min       Q1   median       Q3   max     mean       sd  n missing
#    12.92136 82.16981 93.42432 98.69432 100 87.95012 16.91576 37       0

perc_time_home2_season <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_home_dry$perc_time_home, Nchelenge5_clusts_shorter2_home_rainy$perc_time_home),
                                              c(rep("dry",length(Nchelenge5_clusts_shorter2_home_dry$perc_time_home)),
                                                rep("rainy",length(Nchelenge5_clusts_shorter2_home_rainy$perc_time_home)))))
colnames(perc_time_home2_season) <- c("values","type")
perc_time_home2_season$values <- as.numeric(as.character(perc_time_home2_season$values))
perc_time_home2_season$type <- factor(perc_time_home2_season$type)

ggplot(data=perc_time_home2_season, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual("Season",values = c("sienna", "skyblue")) +
  scale_color_manual("Season",values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,101))
##
kruskal.test(perc_time_home2_season$values,perc_time_home2_season$type)
# p =0.12
### Average Direction of Travel #####
Nchelenge5_clusts_shorter2_no_home$rainy <- as.factor(Nchelenge5_clusts_shorter2_no_home$rainy)
temp <- Nchelenge5_clusts_shorter2_no_home[!duplicated(Nchelenge5_clusts_shorter2_no_home$partid),]
temp2 <- temp[which(!(is.na(temp$rainy))),]

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=rainy), position=position_dodge2(preserve = "single"), width=4) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Rainy Season", values=c("sienna","skyblue"))
#

favstats(temp2$avg_dir_from_home_bearing~temp2$rainy)
#   temp2$rainy       min       Q1   median       Q3      max     mean       sd  n missing
# 1           0 0.3598077 139.3450 225.0665 292.5092 345.7862 207.2764 100.15379 36       0
# 2           1 3.8889726 120.3551 180.2082 288.0634 358.9609 191.8799  91.73408 35       0
favstats(temp2$avg_dir_from_home_mag~temp2$rainy)
#   temp2$rainy       min        Q1   median       Q3       max    mean        sd  n missing
# 1           0 0.1034942 0.5723548 1.444297 2.668043 125.22763 7.443284 23.277144 36       0
# 2           1 0.1811803 0.6604286 1.682104 8.037747  21.89417 4.975965  6.274479 35       0
temp3 <- temp2[which(temp2$avg_dir_from_home_mag < 65),] ## 164 of 176 parts

favstats(temp3$avg_dir_from_home_bearing~temp3$rainy)
#   temp3$rainy       min       Q1   median       Q3      max     mean       sd  n missing
# 1           0 0.3598077 140.7411 227.5121 293.9137 345.7862 210.1898 102.27924 34       0
# 2           1 3.8889726 120.3551 180.2082 288.0634 358.9609 191.8799  91.73408 35       0
favstats(temp3$avg_dir_from_home_mag~temp3$rainy)
#   temp3$rainy       min        Q1   median       Q3      max     mean       sd  n missing
# 1           0 0.1034942 0.5585447 1.390443 2.439780 20.49169 2.178735 3.488989 34       0
# 2           1 0.1811803 0.6604286 1.682104 8.037747 21.89417 4.975965 6.274479 35       0
ggplot(data=temp3, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=rainy), position=position_dodge2(preserve = "single"), width=2) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Rainy Season", values=c("sienna","skyblue"))




temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)] <- -(360-temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)])
ggplot(data=temp2, aes(x=avg_dir_from_home_values, y=..count.., fill=rainy)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~rainy, labeller = labeller(rainy = c("0"="Dry", "1"="Rainy")))+
  scale_fill_manual(values=c("sienna","skyblue")) 

#



temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values < 0)] <- (360+temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values < 0)])

## average direction calc
temp2$values2 <- 90 - temp2$avg_dir_from_home_values
temp2$values2[which(temp2$values2 < 0)] <- temp2$values2[which(temp2$values2 < 0)] + 360
temp2$values3 <- deg2rad(temp2$values2)

temp2$avg_value <- as.numeric(NA)
for(ii in 1:nlevels(temp2$rainy)){
  level <- levels(temp2$rainy)[ii]
  temp <- temp2[which(temp2$rainy == level),]
  sum_cos_dist <- sum(cos(temp$values3))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$values3))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  temp2$avg_value[which(temp2$rainy == level)] <- avg_deg
}

temp2$avg_value2 <- 90 - temp2$avg_value
temp2$avg_value2[which(temp2$avg_value2 < 0)] <- temp2$avg_value2[which(temp2$avg_value2 < 0)] + 360
favstats(temp2$avg_value2~temp2$rainy)
## dry: 280.5 --> W
## rainy: 126.1 --> SE
kruskal.test(temp2$avg_value2,temp2$rainy)

### time spent at clusters vs elsewhere (travel, less important places, etc) ####
# percent time in clusters
Nchelenge5_clusts_shorter2_dry <- Nchelenge5_clusts_shorter2[which(Nchelenge5_clusts_shorter2$rainy == 0),]
Nchelenge5_clusts_shorter2_rainy <- Nchelenge5_clusts_shorter2[which(Nchelenge5_clusts_shorter2$rainy == 1),]
Nchelenge5_clusts_shorter2_dry$perc_time_in_all_clusts <- Nchelenge5_clusts_shorter2_dry$part_trips_time_new/Nchelenge5_clusts_shorter2_dry$total_time
Nchelenge5_clusts_shorter2_rainy$perc_time_in_all_clusts <- Nchelenge5_clusts_shorter2_rainy$part_trips_time_new/Nchelenge5_clusts_shorter2_rainy$total_time

summary(Nchelenge5_clusts_shorter2_dry$perc_time_in_all_clusts[!duplicated(Nchelenge5_clusts_shorter2_dry$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6629  0.8500  0.9261  0.8949  0.9534  0.9955 
summary(Nchelenge5_clusts_shorter2_rainy$perc_time_in_all_clusts[!duplicated(Nchelenge5_clusts_shorter2_rainy$partid)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8287  0.8902  0.9373  0.9318  0.9765  0.9950 

# time outside clusters
summary(Nchelenge5_clusts_shorter2_dry$total_time[!duplicated(Nchelenge5_clusts_shorter2_dry$partid)] - Nchelenge5_clusts_shorter2_dry$part_trips_time_new[!duplicated(Nchelenge5_clusts_shorter2_dry$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.1395  1.1869  2.2344  3.2449  4.6457  9.8070 
summary(Nchelenge5_clusts_shorter2_rainy$total_time[!duplicated(Nchelenge5_clusts_shorter2_rainy$partid)] - Nchelenge5_clusts_shorter2_rainy$part_trips_time_new[!duplicated(Nchelenge5_clusts_shorter2_rainy$partid)])
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max. 
# 0.2383  0.8182  1.8176  2.4824  3.9561  6.6361 
#
#
### Time at locations ####
favstats(Nchelenge5_clusts_shorter2_no_home_dry$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2425 1.080556 2.281944 5.254722 630.46 13.47748 51.19154 357       0
favstats(Nchelenge5_clusts_shorter2_no_home_rainy$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2613889 1.000556 2.200278 4.781181 784.0228 15.09108 65.65735 286       0

### PER PART ###
clust_time_part_median_dry <- (favstats(Nchelenge5_clusts_shorter2_no_home_dry$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_dry$partid)$median)*24
clust_time_part_median_rainy <- (favstats(Nchelenge5_clusts_shorter2_no_home_rainy$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_rainy$partid)$median)*24
favstats(clust_time_part_median_dry)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4605556 1.676181 2.179028 3.204653 17.23194 3.143561 3.172096 36       0
favstats(clust_time_part_median_rainy)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.8686111 1.48125 2.425139 3.184792 12.90486 2.790373 2.162296 35       0

clust_time_part_rainy <- as.data.frame(cbind(c(clust_time_part_median_dry, clust_time_part_median_rainy),
                                             c(rep("dry",length(clust_time_part_median_dry)),
                                               rep("rainy",length(clust_time_part_median_rainy)))))
colnames(clust_time_part_rainy) <- c("values","type")
clust_time_part_rainy$values <- as.numeric(as.character(clust_time_part_rainy$values))
clust_time_part_rainy$type <- factor(clust_time_part_rainy$type)
kruskal.test(clust_time_part_rainy$values,clust_time_part_rainy$type)

ggplot(data=clust_time_part_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-1,20)) +
  scale_x_continuous(breaks=c(seq(0,20,2)))
#
#
#
####### Biting time metrics #####
### number locations outside home ####
Nchelenge5_clusts_shorter2_no_home_biting_places <- Nchelenge5_clusts_shorter2_no_home[which(!(is.na(Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new_biting_time))),]
Nchelenge5_clusts_shorter2_no_home_biting_places2 <- Nchelenge5_clusts_shorter2_no_home_biting_places[which((Nchelenge5_clusts_shorter2_no_home_biting_places$clust_trips_time_new_biting_time*24) >= 0.5),]

Nchelenge5_clusts_shorter2_no_home_biting_places2_dry <- Nchelenge5_clusts_shorter2_no_home_biting_places2[which(Nchelenge5_clusts_shorter2_no_home_biting_places2$rainy == 0),]
Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy <- Nchelenge5_clusts_shorter2_no_home_biting_places2[which(Nchelenge5_clusts_shorter2_no_home_biting_places2$rainy == 1),]

loc_counts_dry <- favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_dry$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_dry$partid)$n
loc_counts_rainy <- favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy$partid)$n
only_home_dry <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_dry$partid))
only_home_rainy <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy$partid))
loc_counts_dry2 <- c(loc_counts_dry,rep(0, only_home_dry))
loc_counts_rainy2 <- c(loc_counts_rainy,rep(0, only_home_rainy))

favstats(loc_counts_dry2)
# min Q1 median Q3 max     mean       sd  n missing
#  0  0      2  3  12 2.135 2.637 37       0
favstats(loc_counts_rainy2)
# min Q1 median Q3 max     mean       sd  n missing
#   0  0      1 1.75   8 1.316 1.861 38       0
summary(as.factor(loc_counts_dry2))
# 0  1  2  3  4  5  6  9 12 
# 13  5  7  4  4  1  1  1  1 
summary(as.factor(loc_counts_rainy2))
# 0  1  2  3  4  6  8 
# 15 13  5  1  1  2  1 
## all those with 0 have biting time at home location (not just nowhere)

locs_count2 <- as.data.frame(cbind(c(loc_counts_dry2, loc_counts_rainy2),
                                   c(rep("dry",length(loc_counts_dry2)),
                                     rep("rainy",length(loc_counts_rainy2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) + 
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-0.35,12.35)) +
  scale_x_continuous(breaks=c(seq(0,12,1)))
#
##
### time per location ####
favstats((Nchelenge5_clusts_shorter2_no_home_biting_places2_dry$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5112 1.804  6.919 16.82 193.3 19.18 35.63 79       0
favstats((Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy$clust_trips_time_new_biting_time)*24) # in hours
#     min    Q1 median    Q3   max  mean    sd   n missing
# 0.5667 1.259  6.056 19.86 285.4 22.98 46.88 50       0

clusts_median_Nchelenge_dry <- (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_dry$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_dry$partid)$median)*24
clusts_median_Nchelenge_rainy <- (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy$partid)$median)*24
favstats(clusts_median_Nchelenge_dry)
# min       Q1   median  Q3      max     mean       sd  n missing
# 1.513 3.894  8.529 19 62.28 14.76 15.2 24       0
favstats(clusts_median_Nchelenge_rainy)
# min       Q1   median  Q3      max     mean       sd  n missing
# 0.8406 1.531   8.26 28.97 82.97 18.98 23.74 23       0

clust_time_part_rainy <- as.data.frame(cbind(c(clusts_median_Nchelenge_dry, clusts_median_Nchelenge_rainy),
                                             c(rep("dry",length(clusts_median_Nchelenge_dry)),
                                               rep("rainy",length(clusts_median_Nchelenge_rainy)))))
colnames(clust_time_part_rainy) <- c("values","type")
clust_time_part_rainy$values <- as.numeric(as.character(clust_time_part_rainy$values))
clust_time_part_rainy$type <- factor(clust_time_part_rainy$type)

ggplot(data=clust_time_part_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,100)) +
  scale_x_continuous(breaks=c(seq(0,100,20)))
#
### number trips ####
Nchelenge5_clust_trips_shorter2_no_home_biting_places <- Nchelenge5_clust_trips_shorter2_no_home[which(!(is.na(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new2_biting_time))),]
## only when >30 minutes during biting time
Nchelenge5_clust_trips_shorter2_no_home_biting_places2 <- Nchelenge5_clust_trips_shorter2_no_home_biting_places[which((Nchelenge5_clust_trips_shorter2_no_home_biting_places$trip_time_new2_biting_time*24) >= 0.5),]
## removes 63 trips, 3 participants

Nchelenge5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2 <- numeric(length=nrow(Nchelenge5_clust_trips_shorter2_no_home_biting_places2))
parts <- levels(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  clusts <- levels(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$clust[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$partid == part)]))
  for(jj in 1: length(clusts)){
    clust <- clusts[jj]
    counts <- Nchelenge5_clust_trips_shorter2_no_home_biting_places2$trip_count_new[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$partid == part & Nchelenge5_clust_trips_shorter2_no_home_biting_places2$clust == clust)]
    counts <- as.factor(counts)
    levels(counts) <- c(1:nlevels(counts))
    Nchelenge5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$partid == part & Nchelenge5_clust_trips_shorter2_no_home_biting_places2$clust == clust)] <- as.numeric(as.character(counts))
  }
}

Nchelenge5_clust_trips_shorter2_no_home_biting_places2b <- Nchelenge5_clust_trips_shorter2_no_home_biting_places2[order(Nchelenge5_clust_trips_shorter2_no_home_biting_places2$partid,
                                                                                                                Nchelenge5_clust_trips_shorter2_no_home_biting_places2$clust,
                                                                                                                -Nchelenge5_clust_trips_shorter2_no_home_biting_places2$trip_count_new2),]
Nchelenge5_clust_trips_shorter2_no_home_biting_places3 <- Nchelenge5_clust_trips_shorter2_no_home_biting_places2b[!duplicated(Nchelenge5_clust_trips_shorter2_no_home_biting_places2b[,c(1,15)]),]

Nchelenge5_clust_trips_shorter2_no_home_biting_places3_dry <- Nchelenge5_clust_trips_shorter2_no_home_biting_places3[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$rainy == 0),]
Nchelenge5_clust_trips_shorter2_no_home_biting_places3_rainy <- Nchelenge5_clust_trips_shorter2_no_home_biting_places3[which(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$rainy == 1),]

favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2  13 2.234 2.486 77       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2)
# min Q1 median Q3 max  mean  sd  n missing
#  1  1      1  2   9 1.898 1.674 49       0

summary(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2))
#   1  2  3  4  6  7  9 10 12 13 
# 45 17  5  2  2  2  1  1  1  1 
summary(as.factor(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2))
#  1  2  3  4  5  7  9 
# 32  5  7  1  2  1  1 

### PER PART ###
trip_count_median_dry <- favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_dry$partid)$median
favstats(trip_count_median_dry)
# min Q1 median Q3 max  mean  sd  n missing
#  1  1    1.5  2 4.5 1.729 0.9998 24       0
trip_count_median_rainy <- favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_rainy$partid)$median
favstats(trip_count_median_rainy)
# min Q1 median Q3 max  mean  sd  n missing
#   1  1      1  2   7 1.804 1.436 23       0

trips_part2_rainy <- as.data.frame(cbind(c(trip_count_median_dry, trip_count_median_rainy),
                                         c(rep("dry",length(trip_count_median_dry)),
                                           rep("rainy",length(trip_count_median_rainy)))))
colnames(trips_part2_rainy) <- c("values","type")
trips_part2_rainy$values <- as.numeric(as.character(trips_part2_rainy$values))
trips_part2_rainy$type <- factor(trips_part2_rainy$type)

ggplot(data=trips_part2_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,10.25)) + 
  scale_x_continuous(breaks=c(seq(0,10,1)))
#
##### percent biting time spent at home vs elsewhere #####
Nchelenge5_clusts_shorter2_biting_places <- Nchelenge5_clusts_shorter2[which(!(is.na(Nchelenge5_clusts_shorter2$clust_trips_time_new_biting_time))),]
Nchelenge5_clusts_shorter2_biting_places$percent_home_biting_time <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_biting_places$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  percent <- Nchelenge5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Nchelenge5_clusts_shorter2_biting_places$partid == part & Nchelenge5_clusts_shorter2_biting_places$home_clust2 == 1)]/
    sum(Nchelenge5_clusts_shorter2_biting_places$clust_trips_time_new_biting_time[which(Nchelenge5_clusts_shorter2_biting_places$partid == part)])
  if(length(percent) == 0){
    percent <- 0
  }
  Nchelenge5_clusts_shorter2_biting_places$percent_home_biting_time[which(Nchelenge5_clusts_shorter2_biting_places$partid == part)] <- percent*100
}
Nchelenge5_clusts_shorter2_biting_places2 <- Nchelenge5_clusts_shorter2_biting_places[!duplicated(Nchelenge5_clusts_shorter2_biting_places$partid),]
homes <- Nchelenge5_clusts_shorter2_home$partid
Nchelenge5_clusts_shorter2_biting_places3 <- Nchelenge5_clusts_shorter2_biting_places[which(Nchelenge5_clusts_shorter2_biting_places$partid %in% homes),]
Nchelenge5_clusts_shorter2_biting_places3b <- Nchelenge5_clusts_shorter2_biting_places3[!duplicated(Nchelenge5_clusts_shorter2_biting_places3$partid),]

Nchelenge5_clusts_shorter2_biting_places3b_dry <- Nchelenge5_clusts_shorter2_biting_places3b[which(Nchelenge5_clusts_shorter2_biting_places3b$rainy == 0),]
Nchelenge5_clusts_shorter2_biting_places3b_rainy <- Nchelenge5_clusts_shorter2_biting_places3b[which(Nchelenge5_clusts_shorter2_biting_places3b$rainy == 1),]

favstats(Nchelenge5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time)
favstats(Nchelenge5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time[which(Nchelenge5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time != 0)])
favstats(Nchelenge5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time)
favstats(Nchelenge5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time[which(Nchelenge5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time != 0)])


perc_time_home2_rainy <- as.data.frame(cbind(c(Nchelenge5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time, Nchelenge5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time),
                                             c(rep("dry",length(Nchelenge5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time)),
                                               rep("rainy",length(Nchelenge5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time)))))
colnames(perc_time_home2_rainy) <- c("values","type")
perc_time_home2_rainy$values <- as.numeric(as.character(perc_time_home2_rainy$values))
perc_time_home2_rainy$type <- factor(perc_time_home2_rainy$type)

ggplot(data=perc_time_home2_rainy, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-1,101))
##

############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
###### ###### TIME OF DAY IN LOCATIONS ###### ######
############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ 
###############
###############
##### CHOMA ####
# -Number of locations where at least 1/4 time spent there during “time of day”
# -Median distance of locations where at least 1/4 time spent there during “time of day”
# -Median number of trips where at least 1/4 of the trip time was during “time of day”
# -Median time per trip for trips where at least 1/4 trip time was during “time of day”
# -Percent time at home (vs elsewhere) during “time of day”
Choma5_clusts_shorter2_morning <- Choma5_clusts_shorter2[which((Choma5_clusts_shorter2$clust_trips_time_new_morning/Choma5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 292 rows
Choma5_clusts_shorter2_midday <- Choma5_clusts_shorter2[which((Choma5_clusts_shorter2$clust_trips_time_new_midday/Choma5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 469 rows
Choma5_clusts_shorter2_evening <- Choma5_clusts_shorter2[which((Choma5_clusts_shorter2$clust_trips_time_new_evening/Choma5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 363 rows
Choma5_clusts_shorter2_night <- Choma5_clusts_shorter2[which((Choma5_clusts_shorter2$clust_trips_time_new_night/Choma5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 37 rows

Choma5_clusts_shorter2_no_home_morning <- Choma5_clusts_shorter2_no_home[which((Choma5_clusts_shorter2_no_home$clust_trips_time_new_morning/Choma5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 240 rows
Choma5_clusts_shorter2_no_home_midday <- Choma5_clusts_shorter2_no_home[which((Choma5_clusts_shorter2_no_home$clust_trips_time_new_midday/Choma5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 448 rows
Choma5_clusts_shorter2_no_home_evening <- Choma5_clusts_shorter2_no_home[which((Choma5_clusts_shorter2_no_home$clust_trips_time_new_evening/Choma5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 317 rows
Choma5_clusts_shorter2_no_home_night <- Choma5_clusts_shorter2_no_home[which((Choma5_clusts_shorter2_no_home$clust_trips_time_new_night/Choma5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 34 rows

Choma5_clust_trips_shorter2_morning <- Choma5_clust_trips_shorter2[which((Choma5_clust_trips_shorter2$trip_time_new2_morning/Choma5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 1210 rows
Choma5_clust_trips_shorter2_midday <- Choma5_clust_trips_shorter2[which((Choma5_clust_trips_shorter2$trip_time_new2_midday/Choma5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 1498 rows
Choma5_clust_trips_shorter2_evening <- Choma5_clust_trips_shorter2[which((Choma5_clust_trips_shorter2$trip_time_new2_evening/Choma5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 1369 rows
Choma5_clust_trips_shorter2_night <- Choma5_clust_trips_shorter2[which((Choma5_clust_trips_shorter2$trip_time_new2_night/Choma5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 305 rows

Choma5_clust_trips_shorter2_no_home_morning <- Choma5_clust_trips_shorter2_no_home[which((Choma5_clust_trips_shorter2_no_home$trip_time_new2_morning/Choma5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 737 rows
Choma5_clust_trips_shorter2_no_home_midday <- Choma5_clust_trips_shorter2_no_home[which((Choma5_clust_trips_shorter2_no_home$trip_time_new2_midday/Choma5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 1119 rows
Choma5_clust_trips_shorter2_no_home_evening <- Choma5_clust_trips_shorter2_no_home[which((Choma5_clust_trips_shorter2_no_home$trip_time_new2_evening/Choma5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 925 rows
Choma5_clust_trips_shorter2_no_home_night <- Choma5_clust_trips_shorter2_no_home[which((Choma5_clust_trips_shorter2_no_home$trip_time_new2_night/Choma5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 186 rows

Choma5_clust_trips_shorter2_no_home3_morning <- Choma5_clust_trips_shorter2_no_home3[which((Choma5_clust_trips_shorter2_no_home3$trip_time_new2_morning/Choma5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 196 rows
Choma5_clust_trips_shorter2_no_home3_midday <- Choma5_clust_trips_shorter2_no_home3[which((Choma5_clust_trips_shorter2_no_home3$trip_time_new2_midday/Choma5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 396 rows
Choma5_clust_trips_shorter2_no_home3_evening <- Choma5_clust_trips_shorter2_no_home3[which((Choma5_clust_trips_shorter2_no_home3$trip_time_new2_evening/Choma5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 278 rows
Choma5_clust_trips_shorter2_no_home3_night <- Choma5_clust_trips_shorter2_no_home3[which((Choma5_clust_trips_shorter2_no_home3$trip_time_new2_night/Choma5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 27 rows
##
###### Split by Time of Day ######
### number of locations (without home) ####
loc_counts_morning <- favstats(Choma5_clusts_shorter2_no_home_morning$clust~Choma5_clusts_shorter2_no_home_morning$partid)$n
nowhere_morning <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_morning$partid))
loc_counts_morning2 <- c(loc_counts_morning,rep(0, nowhere_morning))
favstats(loc_counts_morning2)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0  2      3  5  11 3.806452 2.798353 62       0
loc_counts_midday <- favstats(Choma5_clusts_shorter2_no_home_midday$clust~Choma5_clusts_shorter2_no_home_midday$partid)$n
nowhere_midday <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_midday$partid))
loc_counts_midday2 <- c(loc_counts_midday,rep(0, nowhere_midday))
favstats(loc_counts_midday2)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0  4      8 10  16 7.241935 3.79251 62       0
loc_counts_evening <- favstats(Choma5_clusts_shorter2_no_home_evening$clust~Choma5_clusts_shorter2_no_home_evening$partid)$n
nowhere_evening <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_evening$partid))
loc_counts_evening2 <- c(loc_counts_evening,rep(0, nowhere_evening))
favstats(loc_counts_evening2)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0  2      5  7  12 5.080645 3.250717 62       0
loc_counts_night <- favstats(Choma5_clusts_shorter2_no_home_night$clust~Choma5_clusts_shorter2_no_home_night$partid)$n
nowhere_night <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_night$partid))
loc_counts_night2 <- c(loc_counts_night,rep(0, nowhere_night))
favstats(loc_counts_night2)
# min   Q1 median     Q3  max     mean        sd  n missing
#    0  0      0  1   3 0.6290323 0.8343815 62       0

timecols <- c("#f8f871","#35c7fc", "#fC6a35", "#6a35fc")
locs_count2 <- as.data.frame(cbind(c(loc_counts_morning2, loc_counts_midday2, loc_counts_evening2, loc_counts_night2),
                                   c(rep("morning",length(loc_counts_morning2)),
                                     rep("midday",length(loc_counts_midday2)),
                                     rep("evening",length(loc_counts_evening2)),
                                     rep("night",length(loc_counts_night2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type, levels=c("morning","midday","evening","night"))
kruskal.test(locs_count2$values, locs_count2$type)
pairwise.wilcox.test(locs_count2$values, locs_count2$type, p.adjust.method = "holm")

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.7) + 
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-0.5,25)) +
  scale_x_continuous(breaks=c(seq(0,25,5)))
#
### distance of locations (without home) ####
### TOTAL ###
favstats(Choma5_clusts_shorter2_no_home_morning$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#0.0343186 1.172297 9.84685 16.03122 111.673 11.43895 14.2462 236       0
favstats(Choma5_clusts_shorter2_no_home_midday$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 0.06074684 1.102274 7.407797 15.9999 79.10662 10.60243 12.67453 449       0
favstats(Choma5_clusts_shorter2_no_home_evening$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
# 0.03054886 1.19378 7.727938 14.96106 112.2062 10.55363 13.68201 315       0
favstats(Choma5_clusts_shorter2_no_home_night$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
#  0.427719 3.600859 11.94389 24.45889 111.673 18.52327 23.17064 39       0

### PER PART ###
medians_km_morning <- (favstats(Choma5_clusts_shorter2_no_home_morning$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_morning$partid)$median)
medians_km_midday <- (favstats(Choma5_clusts_shorter2_no_home_midday$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_midday$partid)$median)
medians_km_evening <- (favstats(Choma5_clusts_shorter2_no_home_evening$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_evening$partid)$median)
medians_km_night <- (favstats(Choma5_clusts_shorter2_no_home_night$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_night$partid)$median)

favstats(medians_km_morning) # km
#     min     Q1  median   Q3   max   mean   sd  n missing
#  0.1507284 1.085636 5.551468 14.9576 43.57393 8.994129 9.785163 59       0
favstats(medians_km_midday) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.2040669 1.067656 5.609259 14.85182 55.97212 9.271377 10.32186 60       0
favstats(medians_km_evening) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.3116685 1.091228 5.693963 13.2815 39.1512 7.796477 7.961501 58       0
favstats(medians_km_night) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
# 0.427719 4.383652 11.07769 25.23175 74.90588 18.38595 19.34751 27       0

dists_km2_time <- as.data.frame(cbind(c(medians_km_morning, medians_km_midday, medians_km_evening, medians_km_night),
                                     c(rep("morning",length(medians_km_morning)),
                                       rep("midday",length(medians_km_midday)),
                                       rep("evening",length(medians_km_evening)),
                                       rep("night",length(medians_km_night)))))
colnames(dists_km2_time) <- c("values","type")
dists_km2_time$values <- as.numeric(as.character(dists_km2_time$values))
dists_km2_time$type <- factor(dists_km2_time$type, levels=c("morning","midday","evening","night"))


shapiro.test(dists_km2_time$values) ## not normal
kruskal.test(dists_km2_time$values, dists_km2_time$type) # p = 0.02
pairwise.wilcox.test(dists_km2_time$values, dists_km2_time$type, p.adjust.method = "holm") ## 
#          morning midday evening
#   midday  1.000   -      -      
#   evening 1.000   1.000  -      
#   night   0.040   0.059  0.021 
## night sig diff from all else

ggplot(data=dists_km2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-2.5,85))
#
###
### number of trips (without home) ####
### TOTAL ###
Choma_number_trips_morning <- Choma5_clust_trips_shorter2_no_home_morning[!duplicated(Choma5_clust_trips_shorter2_no_home_morning[,c(1,15)]),c(1,15)]
Choma_number_trips_morning$trips <- as.numeric(0, length(Choma_number_trips_morning$partid))
for(ii in 1:nlevels(as.factor(Choma_number_trips_morning$partid))){
  part <- levels(as.factor(Choma_number_trips_morning$partid))[ii]
  temp <- Choma5_clust_trips_shorter2_no_home_morning[which(Choma5_clust_trips_shorter2_no_home_morning$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Choma_number_trips_morning$trips[which(Choma_number_trips_morning$partid == part & Choma_number_trips_morning$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

Choma_number_trips_midday <- Choma5_clust_trips_shorter2_no_home_midday[!duplicated(Choma5_clust_trips_shorter2_no_home_midday[,c(1,15)]),c(1,15)]
Choma_number_trips_midday$trips <- as.numeric(0, length(Choma_number_trips_midday$partid))
for(ii in 1:nlevels(as.factor(Choma_number_trips_midday$partid))){
  part <- levels(as.factor(Choma_number_trips_midday$partid))[ii]
  temp <- Choma5_clust_trips_shorter2_no_home_midday[which(Choma5_clust_trips_shorter2_no_home_midday$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Choma_number_trips_midday$trips[which(Choma_number_trips_midday$partid == part & Choma_number_trips_midday$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

Choma_number_trips_evening <- Choma5_clust_trips_shorter2_no_home_evening[!duplicated(Choma5_clust_trips_shorter2_no_home_evening[,c(1,15)]),c(1,15)]
Choma_number_trips_evening$trips <- as.numeric(0, length(Choma_number_trips_evening$partid))
for(ii in 1:nlevels(as.factor(Choma_number_trips_evening$partid))){
  part <- levels(as.factor(Choma_number_trips_evening$partid))[ii]
  temp <- Choma5_clust_trips_shorter2_no_home_evening[which(Choma5_clust_trips_shorter2_no_home_evening$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Choma_number_trips_evening$trips[which(Choma_number_trips_evening$partid == part & Choma_number_trips_evening$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

Choma_number_trips_night <- Choma5_clust_trips_shorter2_no_home_night[!duplicated(Choma5_clust_trips_shorter2_no_home_night[,c(1,15)]),c(1,15)]
Choma_number_trips_night$trips <- as.numeric(0, length(Choma_number_trips_night$partid))
for(ii in 1:nlevels(as.factor(Choma_number_trips_night$partid))){
  part <- levels(as.factor(Choma_number_trips_night$partid))[ii]
  temp <- Choma5_clust_trips_shorter2_no_home_night[which(Choma5_clust_trips_shorter2_no_home_night$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Choma_number_trips_night$trips[which(Choma_number_trips_night$partid == part & Choma_number_trips_night$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

### TOTAL ###
favstats(Choma_number_trips_morning$trips)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      1  2  26 2.503497 3.637763 286       0
favstats(Choma_number_trips_midday$trips)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      1  2  31 2.208417 2.897279 499       0
favstats(Choma_number_trips_evening$trips)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      1  2  27 2.378307 3.431367 378       0
favstats(Choma_number_trips_night$trips)
# min Q1 median  Q3 max     mean       sd     n missing
#    1  1      2  3  16 3.257576 3.659805 66       0

### PER PART ###
trips_part_median_morning <- favstats(Choma_number_trips_morning$trips ~ Choma_number_trips_morning$partid)$median
trips_part_median_midday <- favstats(Choma_number_trips_midday$trips ~ Choma_number_trips_midday$partid)$median
trips_part_median_evening <- favstats(Choma_number_trips_evening$trips ~ Choma_number_trips_evening$partid)$median
trips_part_median_night <- favstats(Choma_number_trips_night$trips ~ Choma_number_trips_night$partid)$median

favstats(trips_part_median_morning)
#   min  Q1   median   Q3      max       mean       sd  n missing
# 1  1      1 1.75   6 1.533898 1.058074 59       0
favstats(trips_part_median_midday)
#   min Q1 median  Q3   max      mean       sd  n missing
#  1  1      1 1.125   4 1.266667 0.5559133 60       0
favstats(trips_part_median_evening)
#   min  Q1   median   Q3      max       mean       sd  n missing
#   1  1      1 1.5   5 1.423729 0.9414477 59       0
favstats(trips_part_median_night)
#   min Q1 median  Q3   max      mean       sd  n missing
#  1  1      2 3.375  16 3.157895 3.104281 38       0

trips_part2_time <- as.data.frame(cbind(c(trips_part_median_morning, trips_part_median_midday, trips_part_median_evening, trips_part_median_night),
                                       c(rep("morning",length(trips_part_median_morning)),
                                         rep("midday",length(trips_part_median_midday)),
                                         rep("evening",length(trips_part_median_evening)),
                                         rep("night",length(trips_part_median_night)))))
colnames(trips_part2_time) <- c("values","type")
trips_part2_time$values <- as.numeric(as.character(trips_part2_time$values))
trips_part2_time$type <- factor(trips_part2_time$type, levels=c("morning","midday","evening","night"))

ggplot(data=trips_part2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=10) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(0.5,17)) + 
  scale_x_continuous(breaks=c(seq(0,20,2)))


## remove night outlier 
trips_part_median_night2 <- trips_part_median_night[which(trips_part_median_night < 15)]
favstats(trips_part_median_night2)
# min Q1 median   Q3 max  mean    sd  n missing
#   1  1      2  3  11 2.705882 2.306505 34       0

trips_part2_time <- as.data.frame(cbind(c(trips_part_median_morning, trips_part_median_midday, trips_part_median_evening, trips_part_median_night2),
                                        c(rep("morning",length(trips_part_median_morning)),
                                          rep("midday",length(trips_part_median_midday)),
                                          rep("evening",length(trips_part_median_evening)),
                                          rep("night",length(trips_part_median_night2)))))
colnames(trips_part2_time) <- c("values","type")
trips_part2_time$values <- as.numeric(as.character(trips_part2_time$values))
trips_part2_time$type <- factor(trips_part2_time$type, levels=c("morning","midday","evening","night"))
kruskal.test(trips_part2_time$values, trips_part2_time$type) # p = 0.02
pairwise.wilcox.test(trips_part2_time$values, trips_part2_time$type, p.adjust.method = "holm") ## 

ggplot(data=trips_part2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=8) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(0,11.5)) + 
  scale_x_continuous(breaks=c(seq(0,14,2)))

#
### time per trip ####
### TOTAL ###
favstats(Choma5_clust_trips_shorter2_no_home_morning$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2402778 0.9072917 2.494583 12.61451 472.56 13.9758 38.25074 716       0
favstats(Choma5_clust_trips_shorter2_no_home_midday$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2402778 0.7869444   1.68 4.298056 472.56 10.23749 38.41552 1102       0
favstats(Choma5_clust_trips_shorter2_no_home_evening$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2461111 0.7198611 1.580556 11.44125 472.56 13.87044 40.03681 899       0
favstats(Choma5_clust_trips_shorter2_no_home_night$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.3363889 10.04292 13.68028 19.20875 341.6667 20.03074 28.41495 215       0

### PER PART ###
trip_time_part_median_morning <- favstats(Choma5_clust_trips_shorter2_no_home_morning$trip_time_new2 ~ Choma5_clust_trips_shorter2_no_home_morning$partid)$median
trip_time_part_median_midday <- favstats(Choma5_clust_trips_shorter2_no_home_midday$trip_time_new2 ~ Choma5_clust_trips_shorter2_no_home_midday$partid)$median
trip_time_part_median_evening <- favstats(Choma5_clust_trips_shorter2_no_home_evening$trip_time_new2 ~ Choma5_clust_trips_shorter2_no_home_evening$partid)$median
trip_time_part_median_night <- favstats(Choma5_clust_trips_shorter2_no_home_night$trip_time_new2 ~ Choma5_clust_trips_shorter2_no_home_night$partid)$median

favstats(trip_time_part_median_morning*24)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.5083333 1.68 3.199861 9.820625 163.1999 10.15451 23.72379 59       0
favstats(trip_time_part_median_midday*24)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.7777778 1.338542 1.782639 2.580799 132.4412 4.61834 17.07932 60       0
favstats(trip_time_part_median_evening*24)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4080556 1.022361 1.441806 3.014653 80.07069 5.191784 11.50374 59       0
favstats(trip_time_part_median_night*24)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.2525 10.4267 13.44014 16.69667 56.34111 14.76761 9.351126 38       0

trip_time_part_median_morning_short <- trip_time_part_median_morning[which(trip_time_part_median_morning*24 < 100)]
favstats(trip_time_part_median_morning_short*24)
# min    Q1 median    Q3   max  mean    sd  n missing
# 0.5083333 1.68 3.160069 8.761285 74.26 7.515797 12.43703 58       0
trip_time_part_median_midday_short <- trip_time_part_median_midday[which(trip_time_part_median_midday*24 < 100)]
favstats(trip_time_part_median_midday_short*24)
# min    Q1 median    Q3   max  mean    sd  n missing
# 0.7777778 1.332361 1.765278 2.550139 25.08014 2.45185 3.202527 59       0
trip_time_part_median_evening_short <- trip_time_part_median_evening[which(trip_time_part_median_evening*24 < 100)] 
favstats(trip_time_part_median_evening_short*24)
# min    Q1 median   Q3   max  mean    sd  n missing
# 0.4080556 1.022361 1.441806 3.014653 80.07069 5.191784 11.50374 59       0
trip_time_part_median_night_short <- trip_time_part_median_night[which(trip_time_part_median_night*24 < 100)] 
favstats(trip_time_part_median_night_short*24)
# min    Q1 median   Q3   max  mean    sd  n missing
# 1.2525 10.4267 13.44014 16.69667 56.34111 14.76761 9.351126 38       0

trip_times_part2_time <- as.data.frame(cbind(c(trip_time_part_median_morning_short*24, trip_time_part_median_midday_short*24, trip_time_part_median_evening_short*24, trip_time_part_median_night_short*24),
                                            c(rep("morning",length(trip_time_part_median_morning_short)),
                                              rep("midday",length(trip_time_part_median_midday_short)),
                                              rep("evening",length(trip_time_part_median_evening_short)),
                                              rep("night",length(trip_time_part_median_night_short)))))
colnames(trip_times_part2_time) <- c("values","type")
trip_times_part2_time$values <- as.numeric(as.character(trip_times_part2_time$values))
trip_times_part2_time$type <- factor(trip_times_part2_time$type, levels=c("morning","midday","evening","night"))
kruskal.test(trip_times_part2_time$values, trip_times_part2_time$type)
pairwise.wilcox.test(trip_times_part2_time$values, trip_times_part2_time$type, p.adjust.method = "holm") ## 


ggplot(data=trip_times_part2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=6) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-2,85)) + 
  scale_x_continuous(breaks=c(seq(0,205,10)))

#
### Percent time at home vs elsewhere ####
# -Percent time at home (vs elsewhere) during “time of day”
perc_home_morning <- (Choma5_clusts_shorter2_home$clust_trips_time_new_morning/Choma5_clusts_shorter2_home$part_trips_time_new_morning)*100
perc_home_midday <- (Choma5_clusts_shorter2_home$clust_trips_time_new_midday/Choma5_clusts_shorter2_home$part_trips_time_new_midday)*100
perc_home_evening <- (Choma5_clusts_shorter2_home$clust_trips_time_new_evening/Choma5_clusts_shorter2_home$part_trips_time_new_evening)*100
perc_home_night <- (Choma5_clusts_shorter2_home$clust_trips_time_new_night/Choma5_clusts_shorter2_home$part_trips_time_new_night)*100
perc_home_morning[which(is.na(perc_home_morning))] <- 0
perc_home_midday[which(is.na(perc_home_midday))] <- 0
perc_home_evening[which(is.na(perc_home_evening))] <- 0
perc_home_night[which(is.na(perc_home_night))] <- 0

favstats(perc_home_morning)
# min       Q1   median       Q3 max    mean       sd  n missing
# 0 23.31684 83.58317 96.44422 100 63.15941 37.66578 59       0
favstats(perc_home_midday)
# min      Q1   median       Q3 max     mean       sd  n missing
# 0.2757676 9.80394 68.42161 90.80981 100 55.62051 38.44583 59       0
favstats(perc_home_evening)
# min       Q1   median       Q3 max     mean       sd  n missing
# 0  12.10349 86.6165 96.06207 100 60.44124 40.39291 59       0
favstats(perc_home_night)
# min       Q1  median  Q3 max     mean       sd  n missing
# 0  8.746301 75.28809 100 100 59.65744 42.41673 59       0

perc_home2_time <- as.data.frame(cbind(c(perc_home_morning, perc_home_midday, perc_home_evening, perc_home_night),
                                        c(rep("morning",length(perc_home_morning)),
                                          rep("midday",length(perc_home_midday)),
                                          rep("evening",length(perc_home_evening)),
                                          rep("night",length(perc_home_night)))))
colnames(perc_home2_time) <- c("values","type")
perc_home2_time$values <- as.numeric(as.character(perc_home2_time$values))
perc_home2_time$type <- factor(perc_home2_time$type, levels=c("morning","midday","evening","night"))

shapiro.test(perc_home2_time$values) ## not normal
kruskal.test(perc_home2_time$values, perc_home2_time$type) # p = 0.29
pairwise.wilcox.test(perc_home2_time$values, perc_home2_time$type, p.adjust.method = "holm") ## 

ggplot(data=perc_home2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-7,107)) + 
  scale_x_continuous(breaks=c(seq(0,100,20)))
##

### Time at locations ####
favstats(Choma5_clusts_shorter2_no_home_morning$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2955556 1.277917 2.916111 11.23556 885.6553 51.59654 151.5331 236       0
favstats(Choma5_clusts_shorter2_no_home_midday$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2686111 1.2 2.399722 5.280278 834.3897 28.19682 118.3521 449       0
favstats(Choma5_clusts_shorter2_no_home_evening$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2461111 1.07625 1.914167 4.955 885.6553 46.88331 155.5203 315       0
favstats(Choma5_clusts_shorter2_no_home_night$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.5927778 4.477778 12.5025 55.52139 885.6553 86.05505 191.109 39       0

### PER PART ###
clust_time_part_median_morning <- (favstats(Choma5_clusts_shorter2_no_home_morning$clust_trips_time_new~Choma5_clusts_shorter2_no_home_morning$partid)$median)*24
clust_time_part_median_midday <- (favstats(Choma5_clusts_shorter2_no_home_midday$clust_trips_time_new~Choma5_clusts_shorter2_no_home_midday$partid)$median)*24
clust_time_part_median_evening <- (favstats(Choma5_clusts_shorter2_no_home_evening$clust_trips_time_new~Choma5_clusts_shorter2_no_home_evening$partid)$median)*24
clust_time_part_median_night <- (favstats(Choma5_clusts_shorter2_no_home_night$clust_trips_time_new~Choma5_clusts_shorter2_no_home_night$partid)$median)*24

favstats(clust_time_part_median_morning)
#       min       Q1   median       Q3      max     mean      sd  n missing
#  0.5083333 1.895139 4.235139 9.717778 373.08 31.90361 84.03038 59       0
favstats(clust_time_part_median_midday)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.019167 1.520556 2.314722 3.468785 374.2801 14.56194 62.66199 60       0
favstats(clust_time_part_median_evening)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.72 1.324931 1.937847 3.681701 418.1556 33.67587 103.2689 58       0
favstats(clust_time_part_median_night)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.5927778 5.908333 23.00306 57.995 444.8483 76.04763 128.9232 27       0

clust_time_part_median_morning2 <- clust_time_part_median_morning[which(clust_time_part_median_morning < 400)]
clust_time_part_median_midday2 <- clust_time_part_median_midday[which(clust_time_part_median_midday < 400)]
clust_time_part_median_evening2 <- clust_time_part_median_evening[which(clust_time_part_median_evening < 400)]
clust_time_part_median_night2 <- clust_time_part_median_night[which(clust_time_part_median_night < 400)]
favstats(clust_time_part_median_morning2)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.5083333 1.895139 4.235139 9.717778 373.08 31.90361 84.03038 59       0
favstats(clust_time_part_median_midday2)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 1.019167 1.520556 2.314722 3.468785 374.2801 14.56194 62.66199 60       0
favstats(clust_time_part_median_evening2)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.72 1.318333 1.927361 3.284722 380.8918 26.93061 90.3818 57       0
favstats(clust_time_part_median_night2)
#       min        Q1   median       Q3      max     mean       sd  n missing
#  0.5927778 5.556806 19.04986 51.24278 355.9693 61.86299 107.8713 26       0


clust_time_part_time <- as.data.frame(cbind(c(clust_time_part_median_morning2, clust_time_part_median_midday2,clust_time_part_median_evening2,clust_time_part_median_night2),
                                           c(rep("morning",length(clust_time_part_median_morning2)),
                                             rep("midday",length(clust_time_part_median_midday2)),
                                             rep("evening",length(clust_time_part_median_evening2)),
                                             rep("night",length(clust_time_part_median_night2)))))
colnames(clust_time_part_time) <- c("values","type")
clust_time_part_time$values <- as.numeric(as.character(clust_time_part_time$values))
clust_time_part_time$type <- factor(clust_time_part_time$type, levels=c("morning","midday","evening","night"))
kruskal.test(clust_time_part_time$values, clust_time_part_time$type) # 
pairwise.wilcox.test(clust_time_part_time$values, clust_time_part_time$type, p.adjust.method = "holm") ## 

ggplot(data=clust_time_part_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 50, alpha=0.8) +
  geom_density(aes(color=type), adjust=20) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-10,400)) +
  scale_x_continuous(breaks=c(seq(0,900,50)))
#
#
#
### Average direction #####
######### morning ###
## average degree and magnitude
Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg <- as.numeric(NA)
Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_no_home_morning$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Choma5_clusts_shorter2_no_home_morning[which(Choma5_clusts_shorter2_no_home_morning$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg[which(Choma5_clusts_shorter2_no_home_morning$partid == part)] <- avg_deg
  Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_mag[which(Choma5_clusts_shorter2_no_home_morning$partid == part)] <- avg_vec_dist/1000
}

Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing <- 90 - Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg
Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing[which(Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing < 0)] <- Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing[which(Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_no_home_morning$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Choma5_clusts_shorter2_no_home_morning[which(Choma5_clusts_shorter2_no_home_morning$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg2[which(Choma5_clusts_shorter2_no_home_morning$partid == part)] <- avg_deg
}

Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_values <- 90 - Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg2
Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_values[which(Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_values < 0)] <- Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_values[which(Choma5_clusts_shorter2_no_home_morning$avg_dir_from_home_values < 0)] + 360


######### midday ###
## average degree and magnitude
Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg <- as.numeric(NA)
Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_no_home_midday$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Choma5_clusts_shorter2_no_home_midday[which(Choma5_clusts_shorter2_no_home_midday$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg[which(Choma5_clusts_shorter2_no_home_midday$partid == part)] <- avg_deg
  Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_mag[which(Choma5_clusts_shorter2_no_home_midday$partid == part)] <- avg_vec_dist/1000
}

Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing <- 90 - Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg
Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing[which(Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing < 0)] <- Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing[which(Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_no_home_midday$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Choma5_clusts_shorter2_no_home_midday[which(Choma5_clusts_shorter2_no_home_midday$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg2[which(Choma5_clusts_shorter2_no_home_midday$partid == part)] <- avg_deg
}

Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_values <- 90 - Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg2
Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_values[which(Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_values < 0)] <- Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_values[which(Choma5_clusts_shorter2_no_home_midday$avg_dir_from_home_values < 0)] + 360


######### evening ###
## average degree and magnitude
Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg <- as.numeric(NA)
Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_no_home_evening$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Choma5_clusts_shorter2_no_home_evening[which(Choma5_clusts_shorter2_no_home_evening$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg[which(Choma5_clusts_shorter2_no_home_evening$partid == part)] <- avg_deg
  Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_mag[which(Choma5_clusts_shorter2_no_home_evening$partid == part)] <- avg_vec_dist/1000
}

Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing <- 90 - Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg
Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing[which(Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing < 0)] <- Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing[which(Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_no_home_evening$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Choma5_clusts_shorter2_no_home_evening[which(Choma5_clusts_shorter2_no_home_evening$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg2[which(Choma5_clusts_shorter2_no_home_evening$partid == part)] <- avg_deg
}

Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_values <- 90 - Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg2
Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_values[which(Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_values < 0)] <- Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_values[which(Choma5_clusts_shorter2_no_home_evening$avg_dir_from_home_values < 0)] + 360


######### night ###
## average degree and magnitude
Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_deg <- as.numeric(NA)
Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_no_home_night$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Choma5_clusts_shorter2_no_home_night[which(Choma5_clusts_shorter2_no_home_night$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_deg[which(Choma5_clusts_shorter2_no_home_night$partid == part)] <- avg_deg
  Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_mag[which(Choma5_clusts_shorter2_no_home_night$partid == part)] <- avg_vec_dist/1000
}

Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing <- 90 - Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_deg
Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing[which(Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing < 0)] <- Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing[which(Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Choma5_clusts_shorter2_no_home_night$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Choma5_clusts_shorter2_no_home_night[which(Choma5_clusts_shorter2_no_home_night$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_deg2[which(Choma5_clusts_shorter2_no_home_night$partid == part)] <- avg_deg
}

Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_values <- 90 - Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_deg2
Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_values[which(Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_values < 0)] <- Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_values[which(Choma5_clusts_shorter2_no_home_night$avg_dir_from_home_values < 0)] + 360


######### all ###
temp_morning <- Choma5_clusts_shorter2_no_home_morning[!duplicated(Choma5_clusts_shorter2_no_home_morning$partid),]
temp_midday <- Choma5_clusts_shorter2_no_home_midday[!duplicated(Choma5_clusts_shorter2_no_home_midday$partid),]
temp_evening <- Choma5_clusts_shorter2_no_home_evening[!duplicated(Choma5_clusts_shorter2_no_home_evening$partid),]
temp_night <- Choma5_clusts_shorter2_no_home_night[!duplicated(Choma5_clusts_shorter2_no_home_night$partid),]

temp_dir <- as.data.frame(cbind(c(temp_morning$avg_dir_from_home_bearing, temp_midday$avg_dir_from_home_bearing,temp_evening$avg_dir_from_home_bearing, temp_night$avg_dir_from_home_bearing),
                                c(temp_morning$avg_dir_from_home_mag, temp_midday$avg_dir_from_home_mag,temp_evening$avg_dir_from_home_mag, temp_night$avg_dir_from_home_mag),
                                c(temp_morning$avg_dir_from_home_values, temp_midday$avg_dir_from_home_values,temp_evening$avg_dir_from_home_values, temp_night$avg_dir_from_home_values),
                                c(rep("morning",nrow(temp_morning)),
                                  rep("midday",nrow(temp_midday)),
                                  rep("evening",nrow(temp_evening)),
                                  rep("night",nrow(temp_night)))))

colnames(temp_dir) <- c("bearing","mag", "values", "type")
temp_dir$bearing <- as.numeric(as.character(temp_dir$bearing))
temp_dir$mag <- as.numeric(as.character(temp_dir$mag))
temp_dir$values <- as.numeric(as.character(temp_dir$values))
temp_dir$type <- factor(temp_dir$type, levels=c("morning","midday","evening","night"))

kruskal.test(temp_dir$values, temp_dir$type) # 
pairwise.wilcox.test(temp_dir$values, temp_dir$type, p.adjust.method = "holm") ## 

ggplot(data=temp_dir, aes(x=bearing, y=mag)) +
  geom_col(aes(fill=type), position=position_dodge2(preserve = "single"), width=1.5) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Time of Day", values=timecols, labels=c("morning", "midday", "evening", "night"))


temp_dir$values[which(temp_dir$values > 348.75)] <- -(360-temp_dir$values[which(temp_dir$values > 348.75)])
ggplot(data=temp_dir, aes(x=values, y=..count.., fill=type)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~type)+
  scale_fill_manual(values=timecols) +
  scale_y_continuous(breaks=c(seq(0,14,3)))
##
###############
###############
##### MUTASA ####
# -Number of locations where at least 1/4 time spent there during “time of day”
# -Median distance of locations where at least 1/4 time spent there during “time of day”
# -Median number of trips where at least 1/4 of the trip time was during “time of day”
# -Median time per trip for trips where at least 1/4 trip time was during “time of day”
# -Percent time at home (vs elsewhere) during “time of day”
Mutasa5_clusts_shorter2_morning <- Mutasa5_clusts_shorter2[which((Mutasa5_clusts_shorter2$clust_trips_time_new_morning/Mutasa5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 292 rows
Mutasa5_clusts_shorter2_midday <- Mutasa5_clusts_shorter2[which((Mutasa5_clusts_shorter2$clust_trips_time_new_midday/Mutasa5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 469 rows
Mutasa5_clusts_shorter2_evening <- Mutasa5_clusts_shorter2[which((Mutasa5_clusts_shorter2$clust_trips_time_new_evening/Mutasa5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 363 rows
Mutasa5_clusts_shorter2_night <- Mutasa5_clusts_shorter2[which((Mutasa5_clusts_shorter2$clust_trips_time_new_night/Mutasa5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 37 rows

Mutasa5_clusts_shorter2_no_home_morning <- Mutasa5_clusts_shorter2_no_home[which((Mutasa5_clusts_shorter2_no_home$clust_trips_time_new_morning/Mutasa5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 240 rows
Mutasa5_clusts_shorter2_no_home_midday <- Mutasa5_clusts_shorter2_no_home[which((Mutasa5_clusts_shorter2_no_home$clust_trips_time_new_midday/Mutasa5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 448 rows
Mutasa5_clusts_shorter2_no_home_evening <- Mutasa5_clusts_shorter2_no_home[which((Mutasa5_clusts_shorter2_no_home$clust_trips_time_new_evening/Mutasa5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 317 rows
Mutasa5_clusts_shorter2_no_home_night <- Mutasa5_clusts_shorter2_no_home[which((Mutasa5_clusts_shorter2_no_home$clust_trips_time_new_night/Mutasa5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 34 rows

Mutasa5_clust_trips_shorter2_morning <- Mutasa5_clust_trips_shorter2[which((Mutasa5_clust_trips_shorter2$trip_time_new2_morning/Mutasa5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 1210 rows
Mutasa5_clust_trips_shorter2_midday <- Mutasa5_clust_trips_shorter2[which((Mutasa5_clust_trips_shorter2$trip_time_new2_midday/Mutasa5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 1498 rows
Mutasa5_clust_trips_shorter2_evening <- Mutasa5_clust_trips_shorter2[which((Mutasa5_clust_trips_shorter2$trip_time_new2_evening/Mutasa5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 1369 rows
Mutasa5_clust_trips_shorter2_night <- Mutasa5_clust_trips_shorter2[which((Mutasa5_clust_trips_shorter2$trip_time_new2_night/Mutasa5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 305 rows

Mutasa5_clust_trips_shorter2_no_home_morning <- Mutasa5_clust_trips_shorter2_no_home[which((Mutasa5_clust_trips_shorter2_no_home$trip_time_new2_morning/Mutasa5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 737 rows
Mutasa5_clust_trips_shorter2_no_home_midday <- Mutasa5_clust_trips_shorter2_no_home[which((Mutasa5_clust_trips_shorter2_no_home$trip_time_new2_midday/Mutasa5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 1119 rows
Mutasa5_clust_trips_shorter2_no_home_evening <- Mutasa5_clust_trips_shorter2_no_home[which((Mutasa5_clust_trips_shorter2_no_home$trip_time_new2_evening/Mutasa5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 925 rows
Mutasa5_clust_trips_shorter2_no_home_night <- Mutasa5_clust_trips_shorter2_no_home[which((Mutasa5_clust_trips_shorter2_no_home$trip_time_new2_night/Mutasa5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 186 rows

Mutasa5_clust_trips_shorter2_no_home3_morning <- Mutasa5_clust_trips_shorter2_no_home3[which((Mutasa5_clust_trips_shorter2_no_home3$trip_time_new2_morning/Mutasa5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 196 rows
Mutasa5_clust_trips_shorter2_no_home3_midday <- Mutasa5_clust_trips_shorter2_no_home3[which((Mutasa5_clust_trips_shorter2_no_home3$trip_time_new2_midday/Mutasa5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 396 rows
Mutasa5_clust_trips_shorter2_no_home3_evening <- Mutasa5_clust_trips_shorter2_no_home3[which((Mutasa5_clust_trips_shorter2_no_home3$trip_time_new2_evening/Mutasa5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 278 rows
Mutasa5_clust_trips_shorter2_no_home3_night <- Mutasa5_clust_trips_shorter2_no_home3[which((Mutasa5_clust_trips_shorter2_no_home3$trip_time_new2_night/Mutasa5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 27 rows
##
###### Split by Time of Day ######
### number of locations (without home) ####
loc_counts_morning <- favstats(Mutasa5_clusts_shorter2_no_home_morning$clust~Mutasa5_clusts_shorter2_no_home_morning$partid)$n
nowhere_morning <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_morning$partid))
loc_counts_morning2 <- c(loc_counts_morning,rep(0, nowhere_morning))
favstats(loc_counts_morning2)
# min   Q1 median     Q3  max     mean        sd  n missing
#   0  1      2  4  20 2.818681 2.748195 182       0
loc_counts_midday <- favstats(Mutasa5_clusts_shorter2_no_home_midday$clust~Mutasa5_clusts_shorter2_no_home_midday$partid)$n
nowhere_midday <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_midday$partid))
loc_counts_midday2 <- c(loc_counts_midday,rep(0, nowhere_midday))
favstats(loc_counts_midday2)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0 2.25      5  8  17 5.774725 4.143666 182       0
loc_counts_evening <- favstats(Mutasa5_clusts_shorter2_no_home_evening$clust~Mutasa5_clusts_shorter2_no_home_evening$partid)$n
nowhere_evening <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_evening$partid))
loc_counts_evening2 <- c(loc_counts_evening,rep(0, nowhere_evening))
favstats(loc_counts_evening2)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0  1      2  3  10 2.302198 2.343887 182       0
loc_counts_night <- favstats(Mutasa5_clusts_shorter2_no_home_night$clust~Mutasa5_clusts_shorter2_no_home_night$partid)$n
nowhere_night <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_night$partid))
loc_counts_night2 <- c(loc_counts_night,rep(0, nowhere_night))
favstats(loc_counts_night2)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0  0      0  1   5 0.5604396 0.9128765 182       0

timecols <- c("#f8f871","#35c7fc", "#fC6a35", "#6a35fc")
locs_count2 <- as.data.frame(cbind(c(loc_counts_morning2, loc_counts_midday2, loc_counts_evening2, loc_counts_night2),
                                   c(rep("morning",length(loc_counts_morning2)),
                                     rep("midday",length(loc_counts_midday2)),
                                     rep("evening",length(loc_counts_evening2)),
                                     rep("night",length(loc_counts_night2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type, levels=c("morning","midday","evening","night"))
kruskal.test(locs_count2$values, locs_count2$type)
pairwise.wilcox.test(locs_count2$values, locs_count2$type, p.adjust.method = "holm")

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) + 
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-0.5,25)) +
  scale_x_continuous(breaks=c(seq(0,25,5)))
#
### distance of locations (without home) ####
### TOTAL ###
favstats(Mutasa5_clusts_shorter2_no_home_morning$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
# 4.190736e-05 0.4847252 1.370478 4.795604 481.4493 20.32668 59.97658 513       0
favstats(Mutasa5_clusts_shorter2_no_home_midday$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 2.745604e-05 0.5065108 1.302102 3.450704 481.4493 12.60466 47.2969 1051       0
favstats(Mutasa5_clusts_shorter2_no_home_evening$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
# 2.745604e-05 0.4341955 1.362307 8.540968 476.716 26.78522 69.64332 419       0
favstats(Mutasa5_clusts_shorter2_no_home_night$hhdist_km_haversine) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 0.02519744 0.8790965 6.800982 44.24605 375.1119 45.69367 80.19448 102       0

### PER PART ###
medians_km_morning <- (favstats(Mutasa5_clusts_shorter2_no_home_morning$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_morning$partid)$median)
medians_km_midday <- (favstats(Mutasa5_clusts_shorter2_no_home_midday$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_midday$partid)$median)
medians_km_evening <- (favstats(Mutasa5_clusts_shorter2_no_home_evening$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_evening$partid)$median)
medians_km_night <- (favstats(Mutasa5_clusts_shorter2_no_home_night$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_night$partid)$median)

favstats(medians_km_morning) # km
#     min     Q1  median   Q3   max   mean   sd  n missing
# 0.04779224 0.711989 1.293131 3.488201 350.0869 17.61082 53.35849 156       0
favstats(medians_km_midday) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.08020229 0.6935236 1.193735 2.2506 382.306 9.164458 38.63548 173       0
favstats(medians_km_evening) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.04074332 0.5251485 1.495456 9.902967 374.371 20.67118 55.50707 142       0
favstats(medians_km_night) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
# 0.02519744 1.113139 5.825538 26.69487 271.4218 36.14896 70.241 66       0

medians_km_morning_short <- medians_km_morning[which(medians_km_morning <250)]
medians_km_midday_short <- medians_km_midday[which(medians_km_midday <250)]
medians_km_evening_short <- medians_km_evening[which(medians_km_evening <250)]
medians_km_night_short <- medians_km_night[which(medians_km_night <250)]
favstats(medians_km_morning_short) # km
#     min     Q1  median   Q3   max   mean   sd  n missing
# 0.04779224 0.6951593 1.250154 3.283573 228.4774 15.46581 46.29461 155       0
favstats(medians_km_midday_short) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.08020229 0.68681 1.178366 2.232653 216.1856 6.99503 26.1239 172       0
favstats(medians_km_evening_short) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.04074332 0.5246714 1.478765 9.70704 240.5051 18.16268 46.93727 141       0
favstats(medians_km_night_short) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
# 0.02519744 1.094882 5.711163 22.74867 225.6351 32.52938 64.28598 65       0


dists_km2_time <- as.data.frame(cbind(c(medians_km_morning_short, medians_km_midday_short, medians_km_evening_short, medians_km_night_short),
                                      c(rep("morning",length(medians_km_morning_short)),
                                        rep("midday",length(medians_km_midday_short)),
                                        rep("evening",length(medians_km_evening_short)),
                                        rep("night",length(medians_km_night_short)))))
colnames(dists_km2_time) <- c("values","type")
dists_km2_time$values <- as.numeric(as.character(dists_km2_time$values))
dists_km2_time$type <- factor(dists_km2_time$type, levels=c("morning","midday","evening","night"))

shapiro.test(dists_km2_time$values) ## not normal
kruskal.test(dists_km2_time$values, dists_km2_time$type) # p < 0.001
pairwise.wilcox.test(dists_km2_time$values, dists_km2_time$type, p.adjust.method = "holm") ## 
#           morning midday  evening
#   midday  0.6391  -       -      
#   evening 0.7319  0.6391  -      
#   night   0.0037  4.3e-05 0.0237 
## night sig diff from all else

ggplot(data=dists_km2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=22) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-2.5,240)) +
  scale_x_continuous(breaks=c(seq(0,250,20)))
#
###
### number of trips (without home) ####
### TOTAL ###
Mutasa_number_trips_morning <- Mutasa5_clust_trips_shorter2_no_home_morning[!duplicated(Mutasa5_clust_trips_shorter2_no_home_morning[,c(1,15)]),c(1,15)]
Mutasa_number_trips_morning$trips <- as.numeric(0, length(Mutasa_number_trips_morning$partid))
for(ii in 1:nlevels(as.factor(Mutasa_number_trips_morning$partid))){
  part <- levels(as.factor(Mutasa_number_trips_morning$partid))[ii]
  temp <- Mutasa5_clust_trips_shorter2_no_home_morning[which(Mutasa5_clust_trips_shorter2_no_home_morning$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Mutasa_number_trips_morning$trips[which(Mutasa_number_trips_morning$partid == part & Mutasa_number_trips_morning$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

Mutasa_number_trips_midday <- Mutasa5_clust_trips_shorter2_no_home_midday[!duplicated(Mutasa5_clust_trips_shorter2_no_home_midday[,c(1,15)]),c(1,15)]
Mutasa_number_trips_midday$trips <- as.numeric(0, length(Mutasa_number_trips_midday$partid))
for(ii in 1:nlevels(as.factor(Mutasa_number_trips_midday$partid))){
  part <- levels(as.factor(Mutasa_number_trips_midday$partid))[ii]
  temp <- Mutasa5_clust_trips_shorter2_no_home_midday[which(Mutasa5_clust_trips_shorter2_no_home_midday$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Mutasa_number_trips_midday$trips[which(Mutasa_number_trips_midday$partid == part & Mutasa_number_trips_midday$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

Mutasa_number_trips_evening <- Mutasa5_clust_trips_shorter2_no_home_evening[!duplicated(Mutasa5_clust_trips_shorter2_no_home_evening[,c(1,15)]),c(1,15)]
Mutasa_number_trips_evening$trips <- as.numeric(0, length(Mutasa_number_trips_evening$partid))
for(ii in 1:nlevels(as.factor(Mutasa_number_trips_evening$partid))){
  part <- levels(as.factor(Mutasa_number_trips_evening$partid))[ii]
  temp <- Mutasa5_clust_trips_shorter2_no_home_evening[which(Mutasa5_clust_trips_shorter2_no_home_evening$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Mutasa_number_trips_evening$trips[which(Mutasa_number_trips_evening$partid == part & Mutasa_number_trips_evening$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

Mutasa_number_trips_night <- Mutasa5_clust_trips_shorter2_no_home_night[!duplicated(Mutasa5_clust_trips_shorter2_no_home_night[,c(1,15)]),c(1,15)]
Mutasa_number_trips_night$trips <- as.numeric(0, length(Mutasa_number_trips_night$partid))
for(ii in 1:nlevels(as.factor(Mutasa_number_trips_night$partid))){
  part <- levels(as.factor(Mutasa_number_trips_night$partid))[ii]
  temp <- Mutasa5_clust_trips_shorter2_no_home_night[which(Mutasa5_clust_trips_shorter2_no_home_night$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Mutasa_number_trips_night$trips[which(Mutasa_number_trips_night$partid == part & Mutasa_number_trips_night$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

### TOTAL ###
favstats(Mutasa_number_trips_morning$trips)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      1  2  17 1.884735 2.072579 642       0
favstats(Mutasa_number_trips_midday$trips)
# min Q1 median  Q3 max     mean       sd     n missing
# 1  1      1  2  40 2.106796 2.7777 1133       0
favstats(Mutasa_number_trips_evening$trips)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      1  2  28 2.019784 2.585051 556       0
favstats(Mutasa_number_trips_night$trips)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      1  1  18 1.932836 2.549356 134       0

### PER PART ###
trips_part_median_morning <- favstats(Mutasa_number_trips_morning$trips ~ Mutasa_number_trips_morning$partid)$median
trips_part_median_midday <- favstats(Mutasa_number_trips_midday$trips ~ Mutasa_number_trips_midday$partid)$median
trips_part_median_evening <- favstats(Mutasa_number_trips_evening$trips ~ Mutasa_number_trips_evening$partid)$median
trips_part_median_night <- favstats(Mutasa_number_trips_night$trips ~ Mutasa_number_trips_night$partid)$median

favstats(trips_part_median_morning)
#   min  Q1   median   Q3      max       mean       sd  n missing
#    1  1      1  2  12 1.754601 1.791068 163       0
favstats(trips_part_median_midday)
#   min Q1 median  Q3   max      mean       sd  n missing
#    1  1      1  2  13 1.637931 1.507572 174       0
favstats(trips_part_median_evening)
#   min  Q1   median   Q3      max       mean       sd  n missing
#    1    1      1  2   8 1.645695 1.237726 151       0
favstats(trips_part_median_night)
#   min Q1 median  Q3   max      mean       sd  n missing
#   1  1      1 1.5  14 1.867089 2.237084 79       0

trips_part2_time <- as.data.frame(cbind(c(trips_part_median_morning, trips_part_median_midday, trips_part_median_evening, trips_part_median_night),
                                        c(rep("morning",length(trips_part_median_morning)),
                                          rep("midday",length(trips_part_median_midday)),
                                          rep("evening",length(trips_part_median_evening)),
                                          rep("night",length(trips_part_median_night)))))
colnames(trips_part2_time) <- c("values","type")
trips_part2_time$values <- as.numeric(as.character(trips_part2_time$values))
trips_part2_time$type <- factor(trips_part2_time$type, levels=c("morning","midday","evening","night"))
kruskal.test(trips_part2_time$values, trips_part2_time$type) # 
pairwise.wilcox.test(trips_part2_time$values, trips_part2_time$type, p.adjust.method = "holm") ## 

ggplot(data=trips_part2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(0,14)) + 
  scale_x_continuous(breaks=c(seq(0,20,2)))
#
### time per trip ####
### TOTAL ###
favstats(Mutasa5_clust_trips_shorter2_no_home_morning$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.48 1.609931 3.405278 7.572986 1218.114 12.18537 57.96391 1210       0
favstats(Mutasa5_clust_trips_shorter2_no_home_midday$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.4413889 1.652778 3.120278 5.523194 572.6725 5.407394 18.00812 2387       0
favstats(Mutasa5_clust_trips_shorter2_no_home_evening$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.48 1.310278 2.434167 7.552639 1218.114 13.37646 61.26869 1123       0
favstats(Mutasa5_clust_trips_shorter2_no_home_night$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.5808333 8.64 11.30167 14.96903 182.8889 13.93304 17.06914 259       0

### PER PART ###
trip_time_part_median_morning <- favstats(Mutasa5_clust_trips_shorter2_no_home_morning$trip_time_new2 ~ Mutasa5_clust_trips_shorter2_no_home_morning$partid)$median
trip_time_part_median_midday <- favstats(Mutasa5_clust_trips_shorter2_no_home_midday$trip_time_new2 ~ Mutasa5_clust_trips_shorter2_no_home_midday$partid)$median
trip_time_part_median_evening <- favstats(Mutasa5_clust_trips_shorter2_no_home_evening$trip_time_new2 ~ Mutasa5_clust_trips_shorter2_no_home_evening$partid)$median
trip_time_part_median_night <- favstats(Mutasa5_clust_trips_shorter2_no_home_night$trip_time_new2 ~ Mutasa5_clust_trips_shorter2_no_home_night$partid)$median

favstats(trip_time_part_median_morning*24)
#       min       Q1   median       Q3      max     mean      sd  n missing
#  0.7536111 2.095278  3.205 5.839306 1218.114 21.70392 107.5261 163       0
favstats(trip_time_part_median_midday*24)
#       min        Q1   median       Q3      max     mean       sd  n missing
#  1.144444 2.28375 3.258264 4.180972 11.89319 3.674966 1.977752 174       0
favstats(trip_time_part_median_evening*24)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.48 1.659444 2.271667 4.439861 1218.114 16.48791 101.044 151       0
favstats(trip_time_part_median_night*24)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.5808333 9.145694  11.76 15.43639 182.8889 15.90908 22.21226 79       0

trip_time_part_median_morning_short <- trip_time_part_median_morning[which(trip_time_part_median_morning*24 < 200)]
favstats(trip_time_part_median_morning_short*24)
# min    Q1 median    Q3   max  mean    sd  n missing
# 0.7536111 2.043194 3.092222 5.496319 123.7092 7.690179 17.67661 159       0
trip_time_part_median_evening_short <- trip_time_part_median_evening[which(trip_time_part_median_evening*24 < 200)] 
favstats(trip_time_part_median_evening_short*24)
# min    Q1 median   Q3   max  mean    sd  n missing
# 0.48 1.654653 2.267708 4.303611 181.4308 8.477067 22.8729 150       0

trip_times_part2_time <- as.data.frame(cbind(c(trip_time_part_median_morning_short*24, trip_time_part_median_midday*24, trip_time_part_median_evening_short*24, trip_time_part_median_night*24),
                                             c(rep("morning",length(trip_time_part_median_morning_short)),
                                               rep("midday",length(trip_time_part_median_midday)),
                                               rep("evening",length(trip_time_part_median_evening_short)),
                                               rep("night",length(trip_time_part_median_night)))))
colnames(trip_times_part2_time) <- c("values","type")
trip_times_part2_time$values <- as.numeric(as.character(trip_times_part2_time$values))
trip_times_part2_time$type <- factor(trip_times_part2_time$type, levels=c("morning","midday","evening","night"))

kruskal.test(trip_times_part2_time$values, trip_times_part2_time$type) # 
pairwise.wilcox.test(trip_times_part2_time$values, trip_times_part2_time$type, p.adjust.method = "holm") ## 


ggplot(data=trip_times_part2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=15) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-2,185)) + 
  scale_x_continuous(breaks=c(seq(0,200,20)))

#
### Percent time at home vs elsewhere ####
# -Percent time at home (vs elsewhere) during “time of day”
perc_home_morning <- (Mutasa5_clusts_shorter2_home$clust_trips_time_new_morning/Mutasa5_clusts_shorter2_home$part_trips_time_new_morning)*100
perc_home_midday <- (Mutasa5_clusts_shorter2_home$clust_trips_time_new_midday/Mutasa5_clusts_shorter2_home$part_trips_time_new_midday)*100
perc_home_evening <- (Mutasa5_clusts_shorter2_home$clust_trips_time_new_evening/Mutasa5_clusts_shorter2_home$part_trips_time_new_evening)*100
perc_home_night <- (Mutasa5_clusts_shorter2_home$clust_trips_time_new_night/Mutasa5_clusts_shorter2_home$part_trips_time_new_night)*100
perc_home_morning[which(is.na(perc_home_morning))] <- 0
perc_home_midday[which(is.na(perc_home_midday))] <- 0
perc_home_evening[which(is.na(perc_home_evening))] <- 0
perc_home_night[which(is.na(perc_home_night))] <- 0

favstats(perc_home_morning)
# min       Q1   median       Q3 max    mean       sd  n missing
# 13.67227 76.70773 91.84939 97.68453 100 84.66177 18.28841 182       0
favstats(perc_home_midday)
# min      Q1   median       Q3 max     mean       sd  n missing
# 1.187644 50.32375 76.16539 90.87911 100 68.50957 26.92507 182       0
favstats(perc_home_evening)
# min       Q1   median       Q3 max     mean       sd  n missing
# 5.072642 78.60017 94.69694 98.98898 100 85.32807 20.30945 182       0
favstats(perc_home_night)
# min       Q1  median  Q3 max     mean       sd  n missing
#  0 77.06598 99.23451 100 100 83.34501 25.98655 182       0

perc_home2_time <- as.data.frame(cbind(c(perc_home_morning, perc_home_midday, perc_home_evening, perc_home_night),
                                       c(rep("morning",length(perc_home_morning)),
                                         rep("midday",length(perc_home_midday)),
                                         rep("evening",length(perc_home_evening)),
                                         rep("night",length(perc_home_night)))))
colnames(perc_home2_time) <- c("values","type")
perc_home2_time$values <- as.numeric(as.character(perc_home2_time$values))
perc_home2_time$type <- factor(perc_home2_time$type, levels=c("morning","midday","evening","night"))

shapiro.test(perc_home2_time$values) ## not normal
kruskal.test(perc_home2_time$values, perc_home2_time$type) # p < 0.001
pairwise.wilcox.test(perc_home2_time$values, perc_home2_time$type, p.adjust.method = "holm") ## 
#            morning midday  evening
#   midday  2.2e-10 -       -      
#   evening 0.09342 7.9e-13 -      
#   night   0.00023 1.3e-13 0.00329

ggplot(data=perc_home2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-7,107)) + 
  scale_x_continuous(breaks=c(seq(0,100,20)))
##

### Time at locations ####
favstats(Mutasa5_clusts_shorter2_no_home_morning$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.4802778 2.811389 5.248889 13.08 1537.291 29.67788 113.1419 513       0
favstats(Mutasa5_clusts_shorter2_no_home_midday$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.5136111 2.481528 4.336667 8.363611 472.5967 12.25376 35.12022 1051       0
favstats(Mutasa5_clusts_shorter2_no_home_evening$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.5191667 2.278889  4.965 13.43444 1537.291 37.78915 130.5769 419       0
favstats(Mutasa5_clusts_shorter2_no_home_night$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
#0.5830556 10.08653 13.73514 23.53937 345.2925 30.70606 51.07284 102       0

### PER PART ###
clust_time_part_median_morning <- (favstats(Mutasa5_clusts_shorter2_no_home_morning$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_morning$partid)$median)*24
clust_time_part_median_midday <- (favstats(Mutasa5_clusts_shorter2_no_home_midday$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_midday$partid)$median)*24
clust_time_part_median_evening <- (favstats(Mutasa5_clusts_shorter2_no_home_evening$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_evening$partid)$median)*24
clust_time_part_median_night <- (favstats(Mutasa5_clusts_shorter2_no_home_night$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_night$partid)$median)*24

favstats(clust_time_part_median_morning)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.7772222 3.415069 6.05875 11.82649 1218.114 39.30192 131.1921 156       0
favstats(clust_time_part_median_midday)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.7202778 3.496944 4.3825 6.480278 110.4386 7.440948 12.29812 173       0
favstats(clust_time_part_median_evening)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.6397222 3.301979 5.541111 12.61528 1218.114 36.0374 124.8831 142       0
favstats(clust_time_part_median_night)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.5830556 10.07215 14.8025 24.47681 182.8889 28.63442 39.05171 66       0

clust_time_part_median_morning2 <- clust_time_part_median_morning[which(clust_time_part_median_morning < 250)]
clust_time_part_median_midday2 <- clust_time_part_median_midday[which(clust_time_part_median_midday < 250)]
clust_time_part_median_evening2 <- clust_time_part_median_evening[which(clust_time_part_median_evening < 250)]
clust_time_part_median_night2 <- clust_time_part_median_night[which(clust_time_part_median_night < 250)]
favstats(clust_time_part_median_morning2)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.7772222 3.404444 5.798333 10.66208 215.0031 18.92336 38.8734 151       0
favstats(clust_time_part_median_midday2)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.7202778 3.496944 4.3825 6.480278 110.4386 7.440948 12.29812 173       0
favstats(clust_time_part_median_evening2)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.6397222 3.163507 5.147639 11.76375 217.7682 18.41618 39.87497 138       0
favstats(clust_time_part_median_night2)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.5830556 10.07215 14.8025 24.47681 182.8889 28.63442 39.05171 66       0

clust_time_part_time <- as.data.frame(cbind(c(clust_time_part_median_morning2, clust_time_part_median_midday2,clust_time_part_median_evening2,clust_time_part_median_night2),
                                            c(rep("morning",length(clust_time_part_median_morning2)),
                                              rep("midday",length(clust_time_part_median_midday2)),
                                              rep("evening",length(clust_time_part_median_evening2)),
                                              rep("night",length(clust_time_part_median_night2)))))
colnames(clust_time_part_time) <- c("values","type")
clust_time_part_time$values <- as.numeric(as.character(clust_time_part_time$values))
clust_time_part_time$type <- factor(clust_time_part_time$type, levels=c("morning","midday","evening","night"))
kruskal.test(clust_time_part_time$values, clust_time_part_time$type) # p < 0.001
pairwise.wilcox.test(clust_time_part_time$values, clust_time_part_time$type, p.adjust.method = "holm") ## 

ggplot(data=clust_time_part_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=20) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-10,225)) +
  scale_x_continuous(breaks=c(seq(0,300,20)))
#
#
#
### Average direction #####
######### morning ###
## average degree and magnitude
Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg <- as.numeric(NA)
Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_no_home_morning$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Mutasa5_clusts_shorter2_no_home_morning[which(Mutasa5_clusts_shorter2_no_home_morning$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg[which(Mutasa5_clusts_shorter2_no_home_morning$partid == part)] <- avg_deg
  Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_mag[which(Mutasa5_clusts_shorter2_no_home_morning$partid == part)] <- avg_vec_dist/1000
}

Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing <- 90 - Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg
Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing[which(Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing < 0)] <- Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing[which(Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_no_home_morning$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Mutasa5_clusts_shorter2_no_home_morning[which(Mutasa5_clusts_shorter2_no_home_morning$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg2[which(Mutasa5_clusts_shorter2_no_home_morning$partid == part)] <- avg_deg
}

Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_values <- 90 - Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg2
Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_values[which(Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_values < 0)] <- Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_values[which(Mutasa5_clusts_shorter2_no_home_morning$avg_dir_from_home_values < 0)] + 360


######### midday ###
## average degree and magnitude
Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg <- as.numeric(NA)
Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_no_home_midday$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Mutasa5_clusts_shorter2_no_home_midday[which(Mutasa5_clusts_shorter2_no_home_midday$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg[which(Mutasa5_clusts_shorter2_no_home_midday$partid == part)] <- avg_deg
  Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_mag[which(Mutasa5_clusts_shorter2_no_home_midday$partid == part)] <- avg_vec_dist/1000
}

Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing <- 90 - Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg
Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing[which(Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing < 0)] <- Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing[which(Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_no_home_midday$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Mutasa5_clusts_shorter2_no_home_midday[which(Mutasa5_clusts_shorter2_no_home_midday$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg2[which(Mutasa5_clusts_shorter2_no_home_midday$partid == part)] <- avg_deg
}

Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_values <- 90 - Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg2
Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_values[which(Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_values < 0)] <- Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_values[which(Mutasa5_clusts_shorter2_no_home_midday$avg_dir_from_home_values < 0)] + 360


######### evening ###
## average degree and magnitude
Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg <- as.numeric(NA)
Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_no_home_evening$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Mutasa5_clusts_shorter2_no_home_evening[which(Mutasa5_clusts_shorter2_no_home_evening$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg[which(Mutasa5_clusts_shorter2_no_home_evening$partid == part)] <- avg_deg
  Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_mag[which(Mutasa5_clusts_shorter2_no_home_evening$partid == part)] <- avg_vec_dist/1000
}

Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing <- 90 - Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg
Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing[which(Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing < 0)] <- Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing[which(Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_no_home_evening$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Mutasa5_clusts_shorter2_no_home_evening[which(Mutasa5_clusts_shorter2_no_home_evening$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg2[which(Mutasa5_clusts_shorter2_no_home_evening$partid == part)] <- avg_deg
}

Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_values <- 90 - Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg2
Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_values[which(Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_values < 0)] <- Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_values[which(Mutasa5_clusts_shorter2_no_home_evening$avg_dir_from_home_values < 0)] + 360


######### night ###
## average degree and magnitude
Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_deg <- as.numeric(NA)
Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_no_home_night$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Mutasa5_clusts_shorter2_no_home_night[which(Mutasa5_clusts_shorter2_no_home_night$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_deg[which(Mutasa5_clusts_shorter2_no_home_night$partid == part)] <- avg_deg
  Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_mag[which(Mutasa5_clusts_shorter2_no_home_night$partid == part)] <- avg_vec_dist/1000
}

Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing <- 90 - Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_deg
Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing[which(Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing < 0)] <- Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing[which(Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Mutasa5_clusts_shorter2_no_home_night$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Mutasa5_clusts_shorter2_no_home_night[which(Mutasa5_clusts_shorter2_no_home_night$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_deg2[which(Mutasa5_clusts_shorter2_no_home_night$partid == part)] <- avg_deg
}

Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_values <- 90 - Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_deg2
Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_values[which(Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_values < 0)] <- Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_values[which(Mutasa5_clusts_shorter2_no_home_night$avg_dir_from_home_values < 0)] + 360


######### all ###
temp_morning <- Mutasa5_clusts_shorter2_no_home_morning[!duplicated(Mutasa5_clusts_shorter2_no_home_morning$partid),]
temp_midday <- Mutasa5_clusts_shorter2_no_home_midday[!duplicated(Mutasa5_clusts_shorter2_no_home_midday$partid),]
temp_evening <- Mutasa5_clusts_shorter2_no_home_evening[!duplicated(Mutasa5_clusts_shorter2_no_home_evening$partid),]
temp_night <- Mutasa5_clusts_shorter2_no_home_night[!duplicated(Mutasa5_clusts_shorter2_no_home_night$partid),]

temp_dir <- as.data.frame(cbind(c(temp_morning$avg_dir_from_home_bearing, temp_midday$avg_dir_from_home_bearing,temp_evening$avg_dir_from_home_bearing, temp_night$avg_dir_from_home_bearing),
                                c(temp_morning$avg_dir_from_home_mag, temp_midday$avg_dir_from_home_mag,temp_evening$avg_dir_from_home_mag, temp_night$avg_dir_from_home_mag),
                                c(temp_morning$avg_dir_from_home_values, temp_midday$avg_dir_from_home_values,temp_evening$avg_dir_from_home_values, temp_night$avg_dir_from_home_values),
                                c(rep("morning",nrow(temp_morning)),
                                  rep("midday",nrow(temp_midday)),
                                  rep("evening",nrow(temp_evening)),
                                  rep("night",nrow(temp_night)))))

colnames(temp_dir) <- c("bearing","mag", "values", "type")
temp_dir$bearing <- as.numeric(as.character(temp_dir$bearing))
temp_dir$mag <- as.numeric(as.character(temp_dir$mag))
temp_dir$values <- as.numeric(as.character(temp_dir$values))
temp_dir$type <- factor(temp_dir$type, levels=c("morning","midday","evening","night"))

kruskal.test(temp_dir$values, temp_dir$type) # p < 0.001
pairwise.wilcox.test(temp_dir$values, temp_dir$type, p.adjust.method = "holm") ## 

ggplot(data=temp_dir, aes(x=bearing, y=mag)) +
  geom_col(aes(fill=type), position=position_dodge2(preserve = "single"), width=1.5) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Time of Day", values=timecols, labels=c("morning", "midday", "evening", "night"))


temp_dir$values[which(temp_dir$values > 348.75)] <- -(360-temp_dir$values[which(temp_dir$values > 348.75)])
ggplot(data=temp_dir, aes(x=values, y=..count.., fill=type)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~type)+
  scale_fill_manual(values=timecols) +
  scale_y_continuous(breaks=c(seq(0,20,5)))
##
###############
###############
##### NCHELENGE ####
# -Number of locations where at least 1/4 time spent there during “time of day”
# -Median distance of locations where at least 1/4 time spent there during “time of day”
# -Median number of trips where at least 1/4 of the trip time was during “time of day”
# -Median time per trip for trips where at least 1/4 trip time was during “time of day”
# -Percent time at home (vs elsewhere) during “time of day”
Nchelenge5_clusts_shorter2_morning <- Nchelenge5_clusts_shorter2[which((Nchelenge5_clusts_shorter2$clust_trips_time_new_morning/Nchelenge5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 292 rows
Nchelenge5_clusts_shorter2_midday <- Nchelenge5_clusts_shorter2[which((Nchelenge5_clusts_shorter2$clust_trips_time_new_midday/Nchelenge5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 469 rows
Nchelenge5_clusts_shorter2_evening <- Nchelenge5_clusts_shorter2[which((Nchelenge5_clusts_shorter2$clust_trips_time_new_evening/Nchelenge5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 363 rows
Nchelenge5_clusts_shorter2_night <- Nchelenge5_clusts_shorter2[which((Nchelenge5_clusts_shorter2$clust_trips_time_new_night/Nchelenge5_clusts_shorter2$clust_trips_time_new) >= 0.25),] # 37 rows

Nchelenge5_clusts_shorter2_no_home_morning <- Nchelenge5_clusts_shorter2_no_home[which((Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new_morning/Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 240 rows
Nchelenge5_clusts_shorter2_no_home_midday <- Nchelenge5_clusts_shorter2_no_home[which((Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new_midday/Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 448 rows
Nchelenge5_clusts_shorter2_no_home_evening <- Nchelenge5_clusts_shorter2_no_home[which((Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new_evening/Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 317 rows
Nchelenge5_clusts_shorter2_no_home_night <- Nchelenge5_clusts_shorter2_no_home[which((Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new_night/Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new) >= 0.25),] # 34 rows

Nchelenge5_clust_trips_shorter2_morning <- Nchelenge5_clust_trips_shorter2[which((Nchelenge5_clust_trips_shorter2$trip_time_new2_morning/Nchelenge5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 1210 rows
Nchelenge5_clust_trips_shorter2_midday <- Nchelenge5_clust_trips_shorter2[which((Nchelenge5_clust_trips_shorter2$trip_time_new2_midday/Nchelenge5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 1498 rows
Nchelenge5_clust_trips_shorter2_evening <- Nchelenge5_clust_trips_shorter2[which((Nchelenge5_clust_trips_shorter2$trip_time_new2_evening/Nchelenge5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 1369 rows
Nchelenge5_clust_trips_shorter2_night <- Nchelenge5_clust_trips_shorter2[which((Nchelenge5_clust_trips_shorter2$trip_time_new2_night/Nchelenge5_clust_trips_shorter2$trip_time_new2) >= 0.25),] # 305 rows

Nchelenge5_clust_trips_shorter2_no_home_morning <- Nchelenge5_clust_trips_shorter2_no_home[which((Nchelenge5_clust_trips_shorter2_no_home$trip_time_new2_morning/Nchelenge5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 737 rows
Nchelenge5_clust_trips_shorter2_no_home_midday <- Nchelenge5_clust_trips_shorter2_no_home[which((Nchelenge5_clust_trips_shorter2_no_home$trip_time_new2_midday/Nchelenge5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 1119 rows
Nchelenge5_clust_trips_shorter2_no_home_evening <- Nchelenge5_clust_trips_shorter2_no_home[which((Nchelenge5_clust_trips_shorter2_no_home$trip_time_new2_evening/Nchelenge5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 925 rows
Nchelenge5_clust_trips_shorter2_no_home_night <- Nchelenge5_clust_trips_shorter2_no_home[which((Nchelenge5_clust_trips_shorter2_no_home$trip_time_new2_night/Nchelenge5_clust_trips_shorter2_no_home$trip_time_new2) >= 0.25),] # 186 rows

Nchelenge5_clust_trips_shorter2_no_home3_morning <- Nchelenge5_clust_trips_shorter2_no_home3[which((Nchelenge5_clust_trips_shorter2_no_home3$trip_time_new2_morning/Nchelenge5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 196 rows
Nchelenge5_clust_trips_shorter2_no_home3_midday <- Nchelenge5_clust_trips_shorter2_no_home3[which((Nchelenge5_clust_trips_shorter2_no_home3$trip_time_new2_midday/Nchelenge5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 396 rows
Nchelenge5_clust_trips_shorter2_no_home3_evening <- Nchelenge5_clust_trips_shorter2_no_home3[which((Nchelenge5_clust_trips_shorter2_no_home3$trip_time_new2_evening/Nchelenge5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 278 rows
Nchelenge5_clust_trips_shorter2_no_home3_night <- Nchelenge5_clust_trips_shorter2_no_home3[which((Nchelenge5_clust_trips_shorter2_no_home3$trip_time_new2_night/Nchelenge5_clust_trips_shorter2_no_home3$trip_time_new2) >= 0.25),] # 27 rows
##
###### Split by Time of Day ######
### number of locations (without home) ####
loc_counts_morning <- favstats(Nchelenge5_clusts_shorter2_no_home_morning$clust~Nchelenge5_clusts_shorter2_no_home_morning$partid)$n
nowhere_morning <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_morning$partid))
loc_counts_morning2 <- c(loc_counts_morning,rep(0, nowhere_morning))
favstats(loc_counts_morning2)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0  1      3  5  10 3.28 2.474519 75       0
loc_counts_midday <- favstats(Nchelenge5_clusts_shorter2_no_home_midday$clust~Nchelenge5_clusts_shorter2_no_home_midday$partid)$n
nowhere_midday <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_midday$partid))
loc_counts_midday2 <- c(loc_counts_midday,rep(0, nowhere_midday))
favstats(loc_counts_midday2)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0  2      4 7.5  14 4.973333 3.723459 75       0
loc_counts_evening <- favstats(Nchelenge5_clusts_shorter2_no_home_evening$clust~Nchelenge5_clusts_shorter2_no_home_evening$partid)$n
nowhere_evening <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_evening$partid))
loc_counts_evening2 <- c(loc_counts_evening,rep(0, nowhere_evening))
favstats(loc_counts_evening2)
# min   Q1 median     Q3  max     mean        sd  n missing
#   0  1      3  5  14 3.56 3.205738 75       0
loc_counts_night <- favstats(Nchelenge5_clusts_shorter2_no_home_night$clust~Nchelenge5_clusts_shorter2_no_home_night$partid)$n
nowhere_night <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_night$partid))
loc_counts_night2 <- c(loc_counts_night,rep(0, nowhere_night))
favstats(loc_counts_night2)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0  0      0  1   5 0.4933333 0.92083 75       0

timecols <- c("#f8f871","#35c7fc", "#fC6a35", "#6a35fc")
locs_count2 <- as.data.frame(cbind(c(loc_counts_morning2, loc_counts_midday2, loc_counts_evening2, loc_counts_night2),
                                   c(rep("morning",length(loc_counts_morning2)),
                                     rep("midday",length(loc_counts_midday2)),
                                     rep("evening",length(loc_counts_evening2)),
                                     rep("night",length(loc_counts_night2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type, levels=c("morning","midday","evening","night"))
kruskal.test(locs_count2$values, locs_count2$type)
pairwise.wilcox.test(locs_count2$values, locs_count2$type, p.adjust.method = "holm")

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) + 
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-0.5,25)) +
  scale_x_continuous(breaks=c(seq(0,25,5)))
#
### distance of locations (without home) ####
### TOTAL ###
Nchelenge5_clusts_shorter2_no_home_morning$hhdist_km_haversine2 <- Nchelenge5_clusts_shorter2_no_home_morning$hhdist_m_haversine2/1000
Nchelenge5_clusts_shorter2_no_home_midday$hhdist_km_haversine2 <- Nchelenge5_clusts_shorter2_no_home_midday$hhdist_m_haversine2/1000
Nchelenge5_clusts_shorter2_no_home_evening$hhdist_km_haversine2 <- Nchelenge5_clusts_shorter2_no_home_evening$hhdist_m_haversine2/1000
Nchelenge5_clusts_shorter2_no_home_night$hhdist_km_haversine2 <- Nchelenge5_clusts_shorter2_no_home_night$hhdist_m_haversine2/1000
favstats(Nchelenge5_clusts_shorter2_no_home_morning$hhdist_km_haversine2) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#  0 0.4305217 1.494248 5.102931 187.2311 8.24913 24.48273 246       0
favstats(Nchelenge5_clusts_shorter2_no_home_midday$hhdist_km_haversine2) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
# 0.009427688 0.3805857 1.498601 6.296255 186.7748 6.759266 17.35868 373       0
favstats(Nchelenge5_clusts_shorter2_no_home_evening$hhdist_km_haversine2) # kilometers
#          min        Q1   median       Q3     max     mean       sd   n missing
#  0 0.3251424 0.7615734 2.549729 187.2311 7.655662 23.51278 267       0
favstats(Nchelenge5_clusts_shorter2_no_home_night$hhdist_km_haversine2) # kilometers
#          min        Q1   median       Q3      max     mean       sd   n missing
#  0.03850304 0.3257667 1.937669 11.49557 146.5504 16.02308 33.95572 37       0

### PER PART ###
medians_km_morning <- (favstats(Nchelenge5_clusts_shorter2_no_home_morning$hhdist_km_haversine2~Nchelenge5_clusts_shorter2_no_home_morning$partid)$median)
medians_km_midday <- (favstats(Nchelenge5_clusts_shorter2_no_home_midday$hhdist_km_haversine2~Nchelenge5_clusts_shorter2_no_home_midday$partid)$median)
medians_km_evening <- (favstats(Nchelenge5_clusts_shorter2_no_home_evening$hhdist_km_haversine2~Nchelenge5_clusts_shorter2_no_home_evening$partid)$median)
medians_km_night <- (favstats(Nchelenge5_clusts_shorter2_no_home_night$hhdist_km_haversine2~Nchelenge5_clusts_shorter2_no_home_night$partid)$median)

favstats(medians_km_morning) # km
#     min     Q1  median   Q3   max   mean   sd  n missing
# 0.1300348 0.7615734 1.554269 2.953418 187.0029 7.219929 25.06183 65       0
favstats(medians_km_midday) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
# 0.06520947 0.5544187 1.372192 4.575711 186.7748 7.001752 24.11 68       0
favstats(medians_km_evening) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
# 0.04575551 0.384108 0.7770706 1.803088 94.58437 5.595076 16.63176 65       0
favstats(medians_km_night) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
# 0.08426891 0.4505873 1.937669 6.917181 88.30963 8.106627 18.55934 23       0

medians_km_morning_short <- medians_km_morning[which(medians_km_morning < 150)]
medians_km_midday_short <- medians_km_midday[which(medians_km_midday < 150)]
favstats(medians_km_morning_short) # km
#     min     Q1  median   Q3   max   mean   sd  n missing
#  0.1300348 0.7383462 1.548137 2.892437 79.8984 4.41082 10.81675 64       0
favstats(medians_km_midday_short) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.06520947 0.5535488 1.361116 4.516621 70.95587 4.318573 9.649831 67       0

dists_km2_time <- as.data.frame(cbind(c(medians_km_morning_short, medians_km_midday_short, medians_km_evening, medians_km_night),
                                      c(rep("morning",length(medians_km_morning_short)),
                                        rep("midday",length(medians_km_midday_short)),
                                        rep("evening",length(medians_km_evening)),
                                        rep("night",length(medians_km_night)))))
colnames(dists_km2_time) <- c("values","type")
dists_km2_time$values <- as.numeric(as.character(dists_km2_time$values))
dists_km2_time$type <- factor(dists_km2_time$type, levels=c("morning","midday","evening","night"))

shapiro.test(dists_km2_time$values) ## not normal
kruskal.test(dists_km2_time$values, dists_km2_time$type) # p < 0.10
pairwise.wilcox.test(dists_km2_time$values, dists_km2_time$type, p.adjust.method = "holm") ## 
#          morning midday evening
#   midday  1.00    -      -      
#   evening 0.13    0.31   -      
#   night   1.00    1.00   0.62   


ggplot(data=dists_km2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=10) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-2.5,96)) +
  scale_x_continuous(breaks=c(seq(0,100,10)))
#
###
### number of trips (without home) ####
### TOTAL ###
Nchelenge_number_trips_morning <- Nchelenge5_clust_trips_shorter2_no_home_morning[!duplicated(Nchelenge5_clust_trips_shorter2_no_home_morning[,c(1,14)]),c(1,14)]
Nchelenge_number_trips_morning$trips <- as.numeric(0, length(Nchelenge_number_trips_morning$partid))
for(ii in 1:nlevels(as.factor(Nchelenge_number_trips_morning$partid))){
  part <- levels(as.factor(Nchelenge_number_trips_morning$partid))[ii]
  temp <- Nchelenge5_clust_trips_shorter2_no_home_morning[which(Nchelenge5_clust_trips_shorter2_no_home_morning$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Nchelenge_number_trips_morning$trips[which(Nchelenge_number_trips_morning$partid == part & Nchelenge_number_trips_morning$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

Nchelenge_number_trips_midday <- Nchelenge5_clust_trips_shorter2_no_home_midday[!duplicated(Nchelenge5_clust_trips_shorter2_no_home_midday[,c(1,14)]),c(1,14)]
Nchelenge_number_trips_midday$trips <- as.numeric(0, length(Nchelenge_number_trips_midday$partid))
for(ii in 1:nlevels(as.factor(Nchelenge_number_trips_midday$partid))){
  part <- levels(as.factor(Nchelenge_number_trips_midday$partid))[ii]
  temp <- Nchelenge5_clust_trips_shorter2_no_home_midday[which(Nchelenge5_clust_trips_shorter2_no_home_midday$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Nchelenge_number_trips_midday$trips[which(Nchelenge_number_trips_midday$partid == part & Nchelenge_number_trips_midday$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

Nchelenge_number_trips_evening <- Nchelenge5_clust_trips_shorter2_no_home_evening[!duplicated(Nchelenge5_clust_trips_shorter2_no_home_evening[,c(1,14)]),c(1,14)]
Nchelenge_number_trips_evening$trips <- as.numeric(0, length(Nchelenge_number_trips_evening$partid))
for(ii in 1:nlevels(as.factor(Nchelenge_number_trips_evening$partid))){
  part <- levels(as.factor(Nchelenge_number_trips_evening$partid))[ii]
  temp <- Nchelenge5_clust_trips_shorter2_no_home_evening[which(Nchelenge5_clust_trips_shorter2_no_home_evening$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Nchelenge_number_trips_evening$trips[which(Nchelenge_number_trips_evening$partid == part & Nchelenge_number_trips_evening$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

Nchelenge_number_trips_night <- Nchelenge5_clust_trips_shorter2_no_home_night[!duplicated(Nchelenge5_clust_trips_shorter2_no_home_night[,c(1,14)]),c(1,14)]
Nchelenge_number_trips_night$trips <- as.numeric(0, length(Nchelenge_number_trips_night$partid))
for(ii in 1:nlevels(as.factor(Nchelenge_number_trips_night$partid))){
  part <- levels(as.factor(Nchelenge_number_trips_night$partid))[ii]
  temp <- Nchelenge5_clust_trips_shorter2_no_home_night[which(Nchelenge5_clust_trips_shorter2_no_home_night$partid == part),]
  for(jj in 1:nlevels(as.factor(temp$clust))){
    clust <- levels(as.factor(temp$clust))[jj]
    Nchelenge_number_trips_night$trips[which(Nchelenge_number_trips_night$partid == part & Nchelenge_number_trips_night$clust == clust)] <- nrow(temp[which(temp$clust == clust),])
  }
}

### TOTAL ###
favstats(Nchelenge_number_trips_morning$trips)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      1  2  13 1.702614 1.486483 306       0
favstats(Nchelenge_number_trips_midday$trips)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      1  2  18 1.736471 1.552907 425       0
favstats(Nchelenge_number_trips_evening$trips)
# min Q1 median  Q3 max     mean       sd     n missing
#   1  1      1  2  17 1.939873 1.882128 316       0
favstats(Nchelenge_number_trips_night$trips)
# min Q1 median  Q3 max     mean       sd     n missing
#  1  1      1  2   8 1.836066 1.603957 61       0

### PER PART ###
trips_part_median_morning <- favstats(Nchelenge_number_trips_morning$trips ~ Nchelenge_number_trips_morning$partid)$median
trips_part_median_midday <- favstats(Nchelenge_number_trips_midday$trips ~ Nchelenge_number_trips_midday$partid)$median
trips_part_median_evening <- favstats(Nchelenge_number_trips_evening$trips ~ Nchelenge_number_trips_evening$partid)$median
trips_part_median_night <- favstats(Nchelenge_number_trips_night$trips ~ Nchelenge_number_trips_night$partid)$median

favstats(trips_part_median_morning)
#   min  Q1   median   Q3      max       mean       sd  n missing
#   1  1      1 1.5   6 1.44697 0.8776226 66       0
favstats(trips_part_median_midday)
#   min Q1 median  Q3   max      mean       sd  n missing
# 1  1      1 1.5 2.5 1.314286 0.4675827 70       0
favstats(trips_part_median_evening)
#   min  Q1   median   Q3      max       mean       sd  n missing
#   1    1      1  2   3 1.393939 0.5651082 66       0
favstats(trips_part_median_night)
#   min Q1 median  Q3   max      mean       sd  n missing
#   1    1      1  2   6 1.833333 1.227089 30       0

trips_part2_time <- as.data.frame(cbind(c(trips_part_median_morning, trips_part_median_midday, trips_part_median_evening, trips_part_median_night),
                                        c(rep("morning",length(trips_part_median_morning)),
                                          rep("midday",length(trips_part_median_midday)),
                                          rep("evening",length(trips_part_median_evening)),
                                          rep("night",length(trips_part_median_night)))))
colnames(trips_part2_time) <- c("values","type")
trips_part2_time$values <- as.numeric(as.character(trips_part2_time$values))
trips_part2_time$type <- factor(trips_part2_time$type, levels=c("morning","midday","evening","night"))
kruskal.test(trips_part2_time$values, trips_part2_time$type) # p < 0.10
pairwise.wilcox.test(trips_part2_time$values, trips_part2_time$type, p.adjust.method = "holm") ## 

ggplot(data=trips_part2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(0,10)) + 
  scale_x_continuous(breaks=c(seq(0,10,1)))
#
### time per trip ####
### TOTAL ###
favstats(Nchelenge5_clust_trips_shorter2_no_home_morning$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2469444 0.7119444 1.979722 4.754167 394.7653 10.81987 36.72924 521       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_midday$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2425 0.666875 1.376944 3.114722 261.8942 5.060933 17.97824 738       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_evening$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2416667 0.5788889 1.176667 2.523611 394.7653 9.496829 36.25693 613       0
favstats(Nchelenge5_clust_trips_shorter2_no_home_night$trip_time_new2*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2963889 9.018889 11.43444 18.56375 244.0533 20.42298 33.19596 112       0

### PER PART ###
trip_time_part_median_morning <- favstats(Nchelenge5_clust_trips_shorter2_no_home_morning$trip_time_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_morning$partid)$median
trip_time_part_median_midday <- favstats(Nchelenge5_clust_trips_shorter2_no_home_midday$trip_time_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_midday$partid)$median
trip_time_part_median_evening <- favstats(Nchelenge5_clust_trips_shorter2_no_home_evening$trip_time_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_evening$partid)$median
trip_time_part_median_night <- favstats(Nchelenge5_clust_trips_shorter2_no_home_night$trip_time_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_night$partid)$median

favstats(trip_time_part_median_morning*24)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.3686111 1.347674 2.164514 4.037847 62.32125 5.131328 9.092199 66       0
favstats(trip_time_part_median_midday*24)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.6044444 1.023472 1.364792 2.534792 11.43167 2.044847 1.792013 70       0
favstats(trip_time_part_median_evening*24)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.3743056 0.8917014 1.483472 2.669826 42.44556 3.593451 6.459548 66       0
favstats(trip_time_part_median_night*24)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.6676389 10.46323 12.38264 17.08226 49.71417 15.60986 10.3214 30       0

trip_time_part_median_morning_short <- trip_time_part_median_morning
trip_time_part_median_evening_short <- trip_time_part_median_evening
trip_time_part_median_night_short <- trip_time_part_median_night

trip_times_part2_time <- as.data.frame(cbind(c(trip_time_part_median_morning*24, trip_time_part_median_midday*24, trip_time_part_median_evening*24, trip_time_part_median_night*24),
                                             c(rep("morning",length(trip_time_part_median_morning)),
                                               rep("midday",length(trip_time_part_median_midday)),
                                               rep("evening",length(trip_time_part_median_evening)),
                                               rep("night",length(trip_time_part_median_night)))))
colnames(trip_times_part2_time) <- c("values","type")
trip_times_part2_time$values <- as.numeric(as.character(trip_times_part2_time$values))
trip_times_part2_time$type <- factor(trip_times_part2_time$type, levels=c("morning","midday","evening","night"))
kruskal.test(trip_times_part2_time$values, trip_times_part2_time$type) # p < 0.10
pairwise.wilcox.test(trip_times_part2_time$values, trip_times_part2_time$type, p.adjust.method = "holm") ## 


ggplot(data=trip_times_part2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=5) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-0.5,61)) + 
  scale_x_continuous(breaks=c(seq(0,105,10)))

#
### Percent time at home vs elsewhere ####
# -Percent time at home (vs elsewhere) during “time of day”
perc_home_morning <- (Nchelenge5_clusts_shorter2_home$clust_trips_time_new_morning/Nchelenge5_clusts_shorter2_home$part_trips_time_new_morning)*100
perc_home_midday <- (Nchelenge5_clusts_shorter2_home$clust_trips_time_new_midday/Nchelenge5_clusts_shorter2_home$part_trips_time_new_midday)*100
perc_home_evening <- (Nchelenge5_clusts_shorter2_home$clust_trips_time_new_evening/Nchelenge5_clusts_shorter2_home$part_trips_time_new_evening)*100
perc_home_night <- (Nchelenge5_clusts_shorter2_home$clust_trips_time_new_night/Nchelenge5_clusts_shorter2_home$part_trips_time_new_night)*100
perc_home_morning[which(is.na(perc_home_morning))] <- 0
perc_home_midday[which(is.na(perc_home_midday))] <- 0
perc_home_evening[which(is.na(perc_home_evening))] <- 0
perc_home_night[which(is.na(perc_home_night))] <- 0

favstats(perc_home_morning)
# min       Q1   median       Q3 max    mean       sd  n missing
# 0 78.61106 91.97487 98.0971 100 83.70568 21.04006 74       0
favstats(perc_home_midday)
# min      Q1   median       Q3 max     mean       sd  n missing
# 3.001261 76.09493 87.04859 97.4376 100 81.35647 20.99456 74       0
favstats(perc_home_evening)
# min       Q1   median       Q3 max     mean       sd  n missing
# 0 79.14579 93.5953 98.45487 100 85.07329 21.83569 74       0
favstats(perc_home_night)
# min       Q1  median  Q3 max     mean       sd  n missing
# 0 78.23391    100 100 100 85.6086 24.29352 74       0

perc_home2_time <- as.data.frame(cbind(c(perc_home_morning, perc_home_midday, perc_home_evening, perc_home_night),
                                       c(rep("morning",length(perc_home_morning)),
                                         rep("midday",length(perc_home_midday)),
                                         rep("evening",length(perc_home_evening)),
                                         rep("night",length(perc_home_night)))))
colnames(perc_home2_time) <- c("values","type")
perc_home2_time$values <- as.numeric(as.character(perc_home2_time$values))
perc_home2_time$type <- factor(perc_home2_time$type, levels=c("morning","midday","evening","night"))

shapiro.test(perc_home2_time$values) ## not normal
kruskal.test(perc_home2_time$values, perc_home2_time$type) # p < 0.001
pairwise.wilcox.test(perc_home2_time$values, perc_home2_time$type, p.adjust.method = "holm") ## 
#            morning midday  evening
#   midday  0.59175 -       -      
#   evening 0.59175 0.23583 -      
#   night   0.01220 0.00091 0.02560
# night sig dif from all else

ggplot(data=perc_home2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-7,107)) + 
  scale_x_continuous(breaks=c(seq(0,100,20)))
##

### Time at locations ####
favstats(Nchelenge5_clusts_shorter2_no_home_morning$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2502778 1.535069 3.287639 9.306875 494.26 20.08173 57.97923 246       0
favstats(Nchelenge5_clusts_shorter2_no_home_midday$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2425 1.138056 2.281944 4.718056 630.46 10.20078 42.63234 373       0
favstats(Nchelenge5_clusts_shorter2_no_home_evening$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 0.2436111 1.166944 2.5075 8.138889 784.0228 24.60696 86.37778 267       0
favstats(Nchelenge5_clusts_shorter2_no_home_night$clust_trips_time_new*24)
#       min        Q1   median       Q3      max     mean       sd    n missing
# 1.241667 11.57889 16.1175 43.51333 244.0533 38.423 52.12267 37       0

### PER PART ###
clust_time_part_median_morning <- (favstats(Nchelenge5_clusts_shorter2_no_home_morning$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_morning$partid)$median)*24
clust_time_part_median_midday <- (favstats(Nchelenge5_clusts_shorter2_no_home_midday$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_midday$partid)$median)*24
clust_time_part_median_evening <- (favstats(Nchelenge5_clusts_shorter2_no_home_evening$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_evening$partid)$median)*24
clust_time_part_median_night <- (favstats(Nchelenge5_clusts_shorter2_no_home_night$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_night$partid)$median)*24

favstats(clust_time_part_median_morning)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.465 2.204583 3.6075 9.560556 195.0947 17.71372 39.28059 65       0
favstats(clust_time_part_median_midday)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.7675 1.574757 2.507986 3.676597 27.53389 4.113985 5.366803 68       0
favstats(clust_time_part_median_evening)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4605556 1.71 2.584167 4.078611 184.3115 13.11953 34.65689 65       0
favstats(clust_time_part_median_night)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 6.103194 12.72764 22.99056 58.68396 244.0533 47.72852 60.86481 23       0

clust_time_part_median_morning2 <- clust_time_part_median_morning[which(clust_time_part_median_morning < 150)]
clust_time_part_median_midday2 <- clust_time_part_median_midday[which(clust_time_part_median_midday < 150)]
clust_time_part_median_evening2 <- clust_time_part_median_evening[which(clust_time_part_median_evening < 150)]
clust_time_part_median_night2 <- clust_time_part_median_night[which(clust_time_part_median_night < 150)]

clust_time_part_time <- as.data.frame(cbind(c(clust_time_part_median_morning2, clust_time_part_median_midday2,clust_time_part_median_evening2,clust_time_part_median_night2),
                                            c(rep("morning",length(clust_time_part_median_morning2)),
                                              rep("midday",length(clust_time_part_median_midday2)),
                                              rep("evening",length(clust_time_part_median_evening2)),
                                              rep("night",length(clust_time_part_median_night2)))))
colnames(clust_time_part_time) <- c("values","type")
clust_time_part_time$values <- as.numeric(as.character(clust_time_part_time$values))
clust_time_part_time$type <- factor(clust_time_part_time$type, levels=c("morning","midday","evening","night"))
favstats(clust_time_part_time$values~clust_time_part_time$type)


kruskal.test(clust_time_part_time$values, clust_time_part_time$type) # p < 0.001
pairwise.wilcox.test(clust_time_part_time$values, clust_time_part_time$type, p.adjust.method = "holm") ## 

ggplot(data=clust_time_part_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  #geom_density(aes(color=type), adjust=30) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-15,275)) +
  scale_x_continuous(breaks=c(seq(0,900,50)))
#
#
#
### Average direction #####
######### morning ###
## average degree and magnitude
Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg <- as.numeric(NA)
Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_no_home_morning$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Nchelenge5_clusts_shorter2_no_home_morning[which(Nchelenge5_clusts_shorter2_no_home_morning$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg[which(Nchelenge5_clusts_shorter2_no_home_morning$partid == part)] <- avg_deg
  Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_mag[which(Nchelenge5_clusts_shorter2_no_home_morning$partid == part)] <- avg_vec_dist/1000
}

Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing <- 90 - Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg
Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing[which(Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing < 0)] <- Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing[which(Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_no_home_morning$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Nchelenge5_clusts_shorter2_no_home_morning[which(Nchelenge5_clusts_shorter2_no_home_morning$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg2[which(Nchelenge5_clusts_shorter2_no_home_morning$partid == part)] <- avg_deg
}

Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_values <- 90 - Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_deg2
Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_values[which(Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_values < 0)] <- Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_values[which(Nchelenge5_clusts_shorter2_no_home_morning$avg_dir_from_home_values < 0)] + 360


######### midday ###
## average degree and magnitude
Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg <- as.numeric(NA)
Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_no_home_midday$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Nchelenge5_clusts_shorter2_no_home_midday[which(Nchelenge5_clusts_shorter2_no_home_midday$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg[which(Nchelenge5_clusts_shorter2_no_home_midday$partid == part)] <- avg_deg
  Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_mag[which(Nchelenge5_clusts_shorter2_no_home_midday$partid == part)] <- avg_vec_dist/1000
}

Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing <- 90 - Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg
Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing[which(Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing < 0)] <- Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing[which(Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_no_home_midday$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Nchelenge5_clusts_shorter2_no_home_midday[which(Nchelenge5_clusts_shorter2_no_home_midday$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg2[which(Nchelenge5_clusts_shorter2_no_home_midday$partid == part)] <- avg_deg
}

Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_values <- 90 - Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_deg2
Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_values[which(Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_values < 0)] <- Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_values[which(Nchelenge5_clusts_shorter2_no_home_midday$avg_dir_from_home_values < 0)] + 360


######### evening ###
## average degree and magnitude
Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg <- as.numeric(NA)
Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_no_home_evening$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Nchelenge5_clusts_shorter2_no_home_evening[which(Nchelenge5_clusts_shorter2_no_home_evening$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg[which(Nchelenge5_clusts_shorter2_no_home_evening$partid == part)] <- avg_deg
  Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_mag[which(Nchelenge5_clusts_shorter2_no_home_evening$partid == part)] <- avg_vec_dist/1000
}

Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing <- 90 - Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg
Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing[which(Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing < 0)] <- Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing[which(Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_no_home_evening$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Nchelenge5_clusts_shorter2_no_home_evening[which(Nchelenge5_clusts_shorter2_no_home_evening$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg2[which(Nchelenge5_clusts_shorter2_no_home_evening$partid == part)] <- avg_deg
}

Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_values <- 90 - Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_deg2
Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_values[which(Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_values < 0)] <- Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_values[which(Nchelenge5_clusts_shorter2_no_home_evening$avg_dir_from_home_values < 0)] + 360


######### night ###
## average degree and magnitude
Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_deg <- as.numeric(NA)
Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_mag <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_no_home_night$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Nchelenge5_clusts_shorter2_no_home_night[which(Nchelenge5_clusts_shorter2_no_home_night$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(temp$hhdist_m_haversine*cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(temp$hhdist_m_haversine*sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_deg[which(Nchelenge5_clusts_shorter2_no_home_night$partid == part)] <- avg_deg
  Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_mag[which(Nchelenge5_clusts_shorter2_no_home_night$partid == part)] <- avg_vec_dist/1000
}

Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing <- 90 - Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_deg
Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing[which(Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing < 0)] <- Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing[which(Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_bearing < 0)] + 360


## average degree only, not accounting for distance from home
Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_deg2 <- as.numeric(NA)
parts <- levels(as.factor(Nchelenge5_clusts_shorter2_no_home_night$partid))
for(ii in 1:length(parts)){
  part <- parts[ii]
  temp <- Nchelenge5_clusts_shorter2_no_home_night[which(Nchelenge5_clusts_shorter2_no_home_night$partid == part),]
  ## dist * cos/sin(theta) gives length of horiz/vert components --> sum & divide by # places 
  ## --> use pythag. to get length of average vector
  ## --> use arctan to get angle of average vector
  sum_cos_dist <- sum(cos(temp$dir_from_home_rad))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$dir_from_home_rad))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_deg2[which(Nchelenge5_clusts_shorter2_no_home_night$partid == part)] <- avg_deg
}

Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_values <- 90 - Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_deg2
Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_values[which(Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_values < 0)] <- Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_values[which(Nchelenge5_clusts_shorter2_no_home_night$avg_dir_from_home_values < 0)] + 360


######### all ###
temp_morning <- Nchelenge5_clusts_shorter2_no_home_morning[!duplicated(Nchelenge5_clusts_shorter2_no_home_morning$partid),]
temp_midday <- Nchelenge5_clusts_shorter2_no_home_midday[!duplicated(Nchelenge5_clusts_shorter2_no_home_midday$partid),]
temp_evening <- Nchelenge5_clusts_shorter2_no_home_evening[!duplicated(Nchelenge5_clusts_shorter2_no_home_evening$partid),]
temp_night <- Nchelenge5_clusts_shorter2_no_home_night[!duplicated(Nchelenge5_clusts_shorter2_no_home_night$partid),]

temp_dir <- as.data.frame(cbind(c(temp_morning$avg_dir_from_home_bearing, temp_midday$avg_dir_from_home_bearing,temp_evening$avg_dir_from_home_bearing, temp_night$avg_dir_from_home_bearing),
                                c(temp_morning$avg_dir_from_home_mag, temp_midday$avg_dir_from_home_mag,temp_evening$avg_dir_from_home_mag, temp_night$avg_dir_from_home_mag),
                                c(temp_morning$avg_dir_from_home_values, temp_midday$avg_dir_from_home_values,temp_evening$avg_dir_from_home_values, temp_night$avg_dir_from_home_values),
                                c(rep("morning",nrow(temp_morning)),
                                  rep("midday",nrow(temp_midday)),
                                  rep("evening",nrow(temp_evening)),
                                  rep("night",nrow(temp_night)))))

colnames(temp_dir) <- c("bearing","mag", "values", "type")
temp_dir$bearing <- as.numeric(as.character(temp_dir$bearing))
temp_dir$mag <- as.numeric(as.character(temp_dir$mag))
temp_dir$values <- as.numeric(as.character(temp_dir$values))
temp_dir$type <- factor(temp_dir$type, levels=c("morning","midday","evening","night"))

kruskal.test(temp_dir$values, temp_dir$type) # p < 0.001
pairwise.wilcox.test(temp_dir$values, temp_dir$type, p.adjust.method = "holm") ## 


ggplot(data=temp_dir, aes(x=bearing, y=mag)) +
  geom_col(aes(fill=type), position=position_dodge2(preserve = "single"), width=2.5) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Time of Day", values=timecols, labels=c("morning", "midday", "evening", "night"))


temp_dir$values[which(temp_dir$values > 348.75)] <- -(360-temp_dir$values[which(temp_dir$values > 348.75)])
ggplot(data=temp_dir, aes(x=values, y=..count.., fill=type)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~type)+
  scale_fill_manual(values=timecols) +
  scale_y_continuous(breaks=c(seq(0,10,2)))
##
###############
###############
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
################### COMBINE ALL 3 ################### 
### number of locations (without home) ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)$n
only_home_Choma <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home$partid))
loc_counts_Choma2 <- c(loc_counts_Choma,rep(0, only_home_Choma))

loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)$n
only_home_Mutasa <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home$partid))
loc_counts_Mutasa2 <- c(loc_counts_Mutasa,rep(0, only_home_Mutasa))

loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)$n
only_home_Nchelenge <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home$partid))
loc_counts_Nchelenge2 <- c(loc_counts_Nchelenge,rep(0, only_home_Nchelenge))

locs_count2 <- c(loc_counts_Choma2, loc_counts_Mutasa2, loc_counts_Nchelenge2)

hist(locs_count2, xlim=c(0,26))
favstats(locs_count2)

#
#### distance for 3 areas ####
medians_km_Nchelenge <- (favstats(Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home$partid)$median)
medians_km_Mutasa <- (favstats(Mutasa5_clusts_shorter2_no_home$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home$partid)$median)
medians_km_Choma <- (favstats(Choma5_clusts_shorter2_no_home$hhdist_km_haversine~Choma5_clusts_shorter2_no_home$partid)$median)

medians_km_Nchelenge2 <- medians_km_Nchelenge[which(medians_km_Nchelenge < 80)]
medians_km_Mutasa2 <- medians_km_Mutasa[which(medians_km_Mutasa < 85)]

favstats(c(medians_km_Nchelenge2, medians_km_Mutasa2, medians_km_Choma))

dists_km2 <- as.data.frame(cbind(c(medians_km_Choma, medians_km_Mutasa2, medians_km_Nchelenge2),
                                 c(rep("median",length(medians_km_Choma)),
                                   rep("median",length(medians_km_Mutasa2)),
                                   rep("median",length(medians_km_Nchelenge2)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 5, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=4) +
  coord_cartesian(xlim=c(0,58))
#
#
#### number of trips for 3 areas ####
trips_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home3$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3$partid)$median
trips_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home3$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3$partid)$median
trips_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home3$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3$partid)$median

favstats(c(trips_part_median_Choma, trips_part_median_Mutasa, trips_part_median_Nchelenge))

trips_part2 <- as.data.frame(cbind(c(trips_part_median_Choma, trips_part_median_Mutasa, trips_part_median_Nchelenge),
                                 c(rep("median",length(trips_part_median_Choma)),
                                   rep("median",length(trips_part_median_Mutasa)),
                                   rep("median",length(trips_part_median_Nchelenge)))))
colnames(trips_part2) <- c("values","type")
trips_part2$values <- as.numeric(as.character(trips_part2$values))
trips_part2$type <- factor(trips_part2$type)

ggplot(data=trips_part2, aes(x=values, y=..density..)) +
  geom_histogram(fill="darkgrey", position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1.5) +
  coord_cartesian(xlim=c(0,13.5))+ 
  scale_x_continuous(breaks=c(seq(0,14,2)))
#
#
#### time per trip for 3 areas ####
trip_time_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home$partid)$median
trip_time_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home$partid)$median
trip_time_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home$partid)$median
## remove outliers ##
trip_time_part_median_Mutasa <- trip_time_part_median_Mutasa[-which(trip_time_part_median_Mutasa > 30)]

favstats(c(trip_time_part_median_Nchelenge,trip_time_part_median_Mutasa,trip_time_part_median_Choma))

trip_times_part2 <- as.data.frame(cbind(c(trip_time_part_median_Choma, trip_time_part_median_Mutasa, trip_time_part_median_Nchelenge),
                                 c(rep("median",length(trip_time_part_median_Choma)),
                                   rep("median",length(trip_time_part_median_Mutasa)),
                                   rep("median",length(trip_time_part_median_Nchelenge)))))
colnames(trip_times_part2) <- c("values","type")
trip_times_part2$values <- as.numeric(as.character(trip_times_part2$values))
trip_times_part2$type <- factor(trip_times_part2$type)

ggplot(data=trip_times_part2, aes(x=values, y=..density..)) +
  geom_histogram(binwidth = 0.5, alpha=0.8, fill="grey") +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(0,12)) +
  scale_x_continuous(breaks=c(seq(0,14,2)))

#
#### Percent time at home ####
favstats(Choma5_clusts_shorter2_home$perc_time_home)
favstats(Mutasa5_clusts_shorter2_home$perc_time_home)
favstats(Nchelenge5_clusts_shorter2_home$perc_time_home)

favstats(c(Choma5_clusts_shorter2_home$perc_time_home, Mutasa5_clusts_shorter2_home$perc_time_home, Nchelenge5_clusts_shorter2_home$perc_time_home))


perc_time_home2 <- as.data.frame(cbind(c(Choma5_clusts_shorter2_home$perc_time_home, Mutasa5_clusts_shorter2_home$perc_time_home, Nchelenge5_clusts_shorter2_home$perc_time_home),
                                       c(rep("median",length(Choma5_clusts_shorter2_home$perc_time_home)),
                                         rep("median",length(Mutasa5_clusts_shorter2_home$perc_time_home)),
                                         rep("median",length(Nchelenge5_clusts_shorter2_home$perc_time_home)))))
colnames(perc_time_home2) <- c("values","type")
perc_time_home2$values <- as.numeric(as.character(perc_time_home2$values))
perc_time_home2$type <- factor(perc_time_home2$type)

ggplot(data=perc_time_home2, aes(x=values, y=..density..)) +
  geom_histogram(fill="grey",binwidth = 10, alpha=0.8) +
  geom_density(color="darkgrey", adjust=1) +
  coord_cartesian(xlim=c(-2,102)) +
  scale_x_continuous(breaks=c(seq(0,100,20)))
#
#### direction ####
temp_Choma <- Choma5_clusts_shorter2_no_home[!duplicated(Choma5_clusts_shorter2_no_home$partid),]
temp_Mutasa <- Mutasa5_clusts_shorter2_no_home[!duplicated(Mutasa5_clusts_shorter2_no_home$partid),]
temp_Mutasa2 <- temp_Mutasa[which(temp_Mutasa$avg_dir_from_home_mag < 75),] ## 164 of 176 parts
temp_Nchelenge <- Nchelenge5_clusts_shorter2_no_home[!duplicated(Nchelenge5_clusts_shorter2_no_home$partid),]
temp_Nchelenge2 <- temp_Nchelenge[which(temp_Nchelenge$avg_dir_from_home_mag < 70),] ## 69 of 71 parts

favstats(c(temp_Mutasa2$avg_dir_from_home_bearing,temp_Choma$avg_dir_from_home_bearing,temp_Nchelenge2$avg_dir_from_home_bearing))

temp_direction <- as.data.frame(cbind(c(temp_Choma$avg_dir_from_home_bearing,
                                        temp_Mutasa2$avg_dir_from_home_bearing,
                                        temp_Nchelenge2$avg_dir_from_home_bearing),
                                      c(temp_Choma$avg_dir_from_home_mag,
                                        temp_Mutasa2$avg_dir_from_home_mag,
                                        temp_Nchelenge2$avg_dir_from_home_mag),
                                      c(rep("Choma",length(temp_Choma$avg_dir_from_home_bearing)),
                                        rep("Mutasa",length(temp_Mutasa2$avg_dir_from_home_bearing)),
                                        rep("Nchelenge",length(temp_Nchelenge2$avg_dir_from_home_bearing)),
                                        rep("Choma",length(temp_Choma$avg_dir_from_home_bearing)),
                                        rep("Mutasa",length(temp_Mutasa2$avg_dir_from_home_bearing)),
                                        rep("Nchelenge",length(temp_Nchelenge2$avg_dir_from_home_bearing)))))
colnames(temp_direction) <- c("bearing","magnitude","type")
temp_direction$bearing <- as.numeric(as.character(temp_direction$bearing))
temp_direction$magnitude <- as.numeric(as.character(temp_direction$magnitude))
temp_direction$type <- factor(temp_direction$type)

ggplot(data=temp_direction, aes(x=bearing, y=magnitude)) +
  geom_col(position=position_dodge2(preserve = "single"), width=3, fill="black") +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() 


ggplot(data=temp_direction, aes(x=bearing, y=..count..)) +
  geom_histogram(position = "dodge", binwidth = 11.25, alpha=0.8) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw()





temp_direction <- as.data.frame(cbind(c(temp_Choma$avg_dir_from_home_values,
                                        temp_Mutasa$avg_dir_from_home_values,
                                        temp_Nchelenge$avg_dir_from_home_values),
                                      c(rep("Choma",length(temp_Choma$avg_dir_from_home_values)),
                                        rep("Mutasa",length(temp_Mutasa$avg_dir_from_home_values)),
                                        rep("Nchelenge",length(temp_Nchelenge$avg_dir_from_home_values)))))
colnames(temp_direction) <- c("values","type")
temp_direction$values <- as.numeric(as.character(temp_direction$values))
temp_direction$type <- factor(temp_direction$type)

temp_direction$values[which(temp_direction$values > 348.75)] <- -(360-temp_direction$values[which(temp_direction$values > 348.75)])
ggplot(data=temp_direction, aes(x=values, y=..count..)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw()



## average direction calc
temp_direction$values2 <- 90 - temp_direction$values
temp_direction$values2[which(temp_direction$values2 < 0)] <- temp_direction$values2[which(temp_direction$values2 < 0)] + 360
temp_direction$values3 <- deg2rad(temp_direction$values2)
sum_cos_dist <- sum(cos(temp_direction$values3))/nrow(temp_direction)
sum_sin_dist <- sum(sin(temp_direction$values3))/nrow(temp_direction)
avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
if(avg_deg < 0){
  avg_deg <- 360+avg_deg
}

avg_deg2 <- 90 - avg_deg
avg_deg2[which(avg_deg2 < 0)] <- avg_deg2[which(avg_deg2 < 0)] + 360
avg_deg2
## 305.01

#
#
#
### Time at locations ####
favstats(Choma5_clusts_shorter2_no_home$clust_trips_time_new*24)
clusts_median_Choma <- favstats((Choma5_clusts_shorter2_no_home$clust_trips_time_new*24)~Choma5_clusts_shorter2_no_home$partid)$median
favstats(clusts_median_Choma)
favstats(Mutasa5_clusts_shorter2_no_home$clust_trips_time_new*24)
clusts_median_Mutasa <- favstats((Mutasa5_clusts_shorter2_no_home$clust_trips_time_new*24)~Mutasa5_clusts_shorter2_no_home$partid)$median
favstats(clusts_median_Mutasa)
clusts_median_Mutasa2 <- clusts_median_Mutasa[which(clusts_median_Mutasa < 150)]
favstats(clusts_median_Mutasa2)
favstats(Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new*24)
clusts_median_Nchelenge <- favstats((Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new*24)~Nchelenge5_clusts_shorter2_no_home$partid)$median
favstats(clusts_median_Nchelenge)
#

favstats(c(clusts_median_Choma,clusts_median_Mutasa2,clusts_median_Nchelenge))
#

clust_time2 <- as.data.frame(cbind(c(clusts_median_Choma, clusts_median_Mutasa2, clusts_median_Nchelenge),
                                   c(rep("Choma",length(clusts_median_Choma)),
                                     rep("Mutasa",length(clusts_median_Mutasa2)),
                                     rep("Nchelenge",length(clusts_median_Nchelenge)))))
colnames(clust_time2) <- c("values","type")
clust_time2$values <- as.numeric(as.character(clust_time2$values))
clust_time2$type <- factor(clust_time2$type)

ggplot(data=clust_time2, aes(x=values, y=..density..)) +
  geom_histogram( position = "dodge", binwidth = 5, alpha=0.8, fill="grey") +
  geom_density(adjust=5, color="darkgrey") + 
  coord_cartesian(xlim=c(-3,113)) +
  scale_x_continuous(breaks=c(seq(0,200,10)))
#### Biting time metrics #####
### number locations ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home_biting_places2$clust~Choma5_clusts_shorter2_no_home_biting_places2$partid)$n
only_home_Choma <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2$partid))
loc_counts_Choma2 <- c(loc_counts_Choma,rep(0, only_home_Choma))

loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home_biting_places2$clust~Mutasa5_clusts_shorter2_no_home_biting_places2$partid)$n
only_home_Mutasa <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2$partid))
loc_counts_Mutasa2 <- c(loc_counts_Mutasa,rep(0, only_home_Mutasa))

loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2$partid)$n
only_home_Nchelenge <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2$partid))
loc_counts_Nchelenge2 <- c(loc_counts_Nchelenge,rep(0, only_home_Nchelenge))

locs_count2 <- as.data.frame(cbind(c(loc_counts_Choma2, loc_counts_Mutasa2, loc_counts_Nchelenge2),
                                   c(rep("Choma",length(loc_counts_Choma2)),
                                     rep("Mutasa",length(loc_counts_Mutasa2)),
                                     rep("Nchelenge",length(loc_counts_Nchelenge2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(position = "dodge", binwidth = 2, alpha=0.8, fill="grey") +
  geom_density(adjust=1.5, col="darkgrey") + 
  coord_cartesian(xlim=c(-0.75,13)) +
  scale_x_continuous(breaks=c(seq(0,14,2)))
##
### time per location ####
clusts_median_Choma <- (favstats(Choma5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2$partid)$median)*24
clusts_median_Mutasa <- (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2$partid)$median)*24
clusts_median_Nchelenge <- (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2$partid)$median)*24
clusts_median_Mutasa2 <- clusts_median_Mutasa[which(clusts_median_Mutasa < 500)]

clust_time2 <- as.data.frame(cbind(c(clusts_median_Choma, clusts_median_Mutasa2, clusts_median_Nchelenge),
                                   c(rep("Choma",length(clusts_median_Choma)),
                                     rep("Mutasa",length(clusts_median_Mutasa2)),
                                     rep("Nchelenge",length(clusts_median_Nchelenge)))))
colnames(clust_time2) <- c("values","type")
clust_time2$values <- as.numeric(as.character(clust_time2$values))
clust_time2$type <- factor(clust_time2$type)

ggplot(data=clust_time2, aes(x=values, y=..density..)) +
  geom_histogram(position = "dodge", binwidth = 20, alpha=0.8, fill="grey") +
  geom_density( adjust=3, col="darkgrey") + 
  coord_cartesian(xlim=c(-18,300)) +
  scale_x_continuous(breaks=c(seq(0,300,50)))
#
### number trips ####
### PER PART ###
trip_count_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3$partid)$median
trip_count_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3$partid)$median
trip_count_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3$partid)$median
trip_count_median_Choma2 <- trip_count_median_Choma[which(trip_count_median_Choma < 20)]

trips_all <- as.data.frame(cbind(c(trip_count_median_Choma2, trip_count_median_Mutasa, trip_count_median_Nchelenge),
                                 c(rep("Choma",length(trip_count_median_Choma2)),
                                   rep("Mutasa",length(trip_count_median_Mutasa)),
                                   rep("Nchelenge",length(trip_count_median_Nchelenge)))))
colnames(trips_all) <- c("values","type")
trips_all$values <- as.numeric(as.character(trips_all$values))
trips_all$type <- factor(trips_all$type)

ggplot(data=trips_all, aes(x=values, y=..density..)) +
  geom_histogram( position = "dodge", binwidth = 1, alpha=0.8, fill="grey") +
  geom_density(adjust=2, col="darkgrey") +
  coord_cartesian(xlim=c(0,13)) +
  scale_x_continuous(breaks=c(seq(0,30,2)))
#

trips_all2 <- as.data.frame(cbind(c(Choma5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2, 
                                    Mutasa5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2, 
                                    Nchelenge5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2),
                                  c(rep("Choma",length(Choma5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2)),
                                    rep("Mutasa",length(Mutasa5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2)),
                                    rep("Nchelenge",length(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2)))))
colnames(trips_all2) <- c("values","type")
trips_all2$values <- as.numeric(as.character(trips_all2$values))
trips_all2$type <- factor(trips_all2$type)

table(as.factor(trips_all2$values))

#
##### percent biting time spent at home vs elsewhere #####
perc_time_home2 <- as.data.frame(cbind(c(Choma5_clusts_shorter2_biting_places3b$percent_home_biting_time, 
                                         Mutasa5_clusts_shorter2_biting_places3b$percent_home_biting_time, 
                                         Nchelenge5_clusts_shorter2_biting_places3b$percent_home_biting_time),
                                       c(rep("Choma",length(Choma5_clusts_shorter2_biting_places3b$percent_home_biting_time)),
                                         rep("Mutasa",length(Mutasa5_clusts_shorter2_biting_places3b$percent_home_biting_time)),
                                         rep("Nchelenge",length(Nchelenge5_clusts_shorter2_biting_places3b$percent_home_biting_time)))))
colnames(perc_time_home2) <- c("values","type")
perc_time_home2$values <- as.numeric(as.character(perc_time_home2$values))
perc_time_home2$type <- factor(perc_time_home2$type)

ggplot(data=perc_time_home2, aes(x=values, y=..density..)) +
  geom_histogram(position = "dodge", binwidth = 10, alpha=0.8, fill="grey") +
  geom_density(adjust=1, col="darkgrey") +
  coord_cartesian(xlim=c(-2,102)) +
  scale_x_continuous(breaks=c(seq(0,100,20)))
##
###############
###############
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
################### COMPARE ACROSS ALL 3 ################### 
### number of locations (without home) ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)$n
only_home_Choma <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home$partid))
loc_counts_Choma2 <- c(loc_counts_Choma,rep(0, only_home_Choma))

loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)$n
only_home_Mutasa <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home$partid))
loc_counts_Mutasa2 <- c(loc_counts_Mutasa,rep(0, only_home_Mutasa))

loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)$n
only_home_Nchelenge <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home$partid))
loc_counts_Nchelenge2 <- c(loc_counts_Nchelenge,rep(0, only_home_Nchelenge))

locs_count2 <- as.data.frame(cbind(c(loc_counts_Choma2, loc_counts_Mutasa2, loc_counts_Nchelenge2),
                                 c(rep("Choma",length(loc_counts_Choma2)),
                                   rep("Mutasa",length(loc_counts_Mutasa2)),
                                   rep("Nchelenge",length(loc_counts_Nchelenge2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)
cols <- brewer.pal(n=4, "Set1")
cols2 <- c(cols[4],SteppedSequential5Steps[c(18,13)])
ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = cols2) +
  scale_color_manual(values = cols2) +
  coord_cartesian(xlim=c(-2,27)) +
  scale_x_continuous(breaks=c(seq(0,26,5)))


shapiro.test(locs_count2$values) ## not normal
kruskal.test(locs_count2$values, locs_count2$type) # p = 1.58e-05
pairwise.wilcox.test(locs_count2$values, locs_count2$type, p.adjust.method = "holm") ## 
#           Choma    Mutasa
# Mutasa    8.8e-06  -     
# Nchelenge 0.0022   0.3297

my_comp <- list(c("Choma", "Mutasa"), c("Choma","Nchelenge"), c("Mutasa","Nchelenge"))
ggplot(data=locs_count2, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = cols2) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 35) 
  
##



descdist(locs_count2$values, discrete = TRUE, boot=1000)
a1 <- fitdist(locs_count2$values, "pois")
a2 <- fitdist(locs_count2$values, "nbinom")
denscomp(list(a1,a2)) ## nbinom better

#
#### distance for 3 areas ####
medians_km_Nchelenge <- (favstats(Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home$partid)$median)
medians_km_Mutasa <- (favstats(Mutasa5_clusts_shorter2_no_home$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home$partid)$median)
medians_km_Choma <- (favstats(Choma5_clusts_shorter2_no_home$hhdist_km_haversine~Choma5_clusts_shorter2_no_home$partid)$median)

medians_km_Nchelenge2 <- medians_km_Nchelenge[which(medians_km_Nchelenge < 80)]
medians_km_Mutasa2 <- medians_km_Mutasa[which(medians_km_Mutasa < 85)]

dists_km2 <- as.data.frame(cbind(c(medians_km_Choma, medians_km_Mutasa2, medians_km_Nchelenge2),
                                 c(rep("Choma",length(medians_km_Choma)),
                                   rep("Mutasa",length(medians_km_Mutasa2)),
                                   rep("Nchelenge",length(medians_km_Nchelenge2)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)

cols <- brewer.pal(n=4, "Set1")
ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=8) + 
  scale_fill_manual(values = cols2) +
  scale_color_manual(values = cols2) +
  coord_cartesian(xlim=c(-2,65)) +
  scale_x_continuous(breaks=c(seq(0,60,20)))

shapiro.test(dists_km2$values) ## not normal
kruskal.test(dists_km2$values, dists_km2$type) # p = 3e-05
pairwise.wilcox.test(dists_km2$values, dists_km2$type, p.adjust.method = "holm") ## 
#           Choma    Mutasa
# Mutasa    5e-05  -     
# Nchelenge  1e-04 0.8   

my_comp <- list(c("Choma", "Mutasa"), c("Choma","Nchelenge"), c("Mutasa","Nchelenge"))
ggplot(data=dists_km2, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = cols2) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.5) +
  stat_compare_means(label.y = 78) 





#######
descdist(dists_km2$values, boot=1000)
a1 <- fitdist(dists_km2$values, "exp")
a2 <- fitdist(dists_km2$values, "gamma")
a3 <- fitdist(dists_km2$values, "weibull")
a4 <- fitdist(dists_km2$values, "lnorm")
a5 <- fitdist(dists_km2$values, "norm")

denscomp(list(a1,a2,a3,a4,a5)) 
ppcomp(list(a1,a2,a3,a4,a5)) ## lnorm better

summary(a4)
#          estimate Std. Error
# meanlog 0.4975619 0.07922799
# sdlog   1.3722691 0.05602252


aa <- glm(values~type, data=dists_km2, family=gaussian(link="log"))



#
#### number of trips for 3 areas ####
trips_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home3$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3$partid)$median
trips_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home3$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3$partid)$median
trips_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home3$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3$partid)$median

dists_km2 <- as.data.frame(cbind(c(trips_part_median_Choma, trips_part_median_Mutasa, trips_part_median_Nchelenge),
                                 c(rep("Choma",length(trips_part_median_Choma)),
                                   rep("Mutasa",length(trips_part_median_Mutasa)),
                                   rep("Nchelenge",length(trips_part_median_Nchelenge)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)
cols <- brewer.pal(n=4, "Set1")

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2.5) +
  scale_fill_manual(values = cols2) +
  scale_color_manual(values = cols2) +
  coord_cartesian(xlim=c(0,13)) +
  scale_x_continuous(breaks=c(seq(0,14,2)))


shapiro.test(dists_km2$values) ## not normal
kruskal.test(dists_km2$values, dists_km2$type) # p = 0.4
pairwise.wilcox.test(dists_km2$values, dists_km2$type, p.adjust.method = "holm") ## 
#           Choma   Mutasa
# Mutasa    0.6     -     
# Nchelenge  0.6    0.8   

my_comp <- list(c("Choma", "Mutasa"), c("Choma","Nchelenge"), c("Mutasa","Nchelenge"))
ggplot(data=dists_km2, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = cols2) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.5) +
  stat_compare_means(label.y = 18) 

#
#### time per trip for 3 areas ####
trip_time_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home$partid)$median
trip_time_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home$partid)$median
trip_time_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home$partid)$median
## remove outliers ##
trip_time_part_median_Mutasa <- trip_time_part_median_Mutasa[-which(trip_time_part_median_Mutasa > 30)]

dists_km2 <- as.data.frame(cbind(c(trip_time_part_median_Choma, trip_time_part_median_Mutasa, trip_time_part_median_Nchelenge),
                                 c(rep("Choma",length(trip_time_part_median_Choma)),
                                   rep("Mutasa",length(trip_time_part_median_Mutasa)),
                                   rep("Nchelenge",length(trip_time_part_median_Nchelenge)))))
colnames(dists_km2) <- c("values","type")
dists_km2$values <- as.numeric(as.character(dists_km2$values))
dists_km2$type <- factor(dists_km2$type)
cols <- brewer.pal(n=4, "Set1")

ggplot(data=dists_km2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=4) +
  scale_fill_manual(values = cols2) +
  scale_color_manual(values = cols2) +
  coord_cartesian(xlim=c(-1,13)) +
  scale_x_continuous(breaks=c(seq(0,15,2)))


shapiro.test(dists_km2$values) ## not normal
kruskal.test(dists_km2$values, dists_km2$type) # p < 2e-16
pairwise.wilcox.test(dists_km2$values, dists_km2$type, p.adjust.method = "holm") ## 
#           Choma   Mutasa
# Mutasa    <2e-16     -     
# Nchelenge  0.3   <2e-16

my_comp <- list(c("Choma", "Mutasa"), c("Choma","Nchelenge"), c("Mutasa","Nchelenge"))
ggplot(data=dists_km2, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = cols2) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.5) +
  stat_compare_means(label.y = 17) 



descdist(dists_km2$values, boot=1000)
a1 <- fitdist(dists_km2$values, "exp")
a2 <- fitdist(dists_km2$values, "gamma")
a3 <- fitdist(dists_km2$values, "weibull")
a4 <- fitdist(dists_km2$values, "lnorm")
a5 <- fitdist(dists_km2$values, "norm")

denscomp(list(a1,a2,a3,a4,a5)) 
ppcomp(list(a1,a2,a3,a4,a5)) ## lnorm better

#
#### Percent time at home ####
favstats(Choma5_clusts_shorter2_home$perc_time_home)
favstats(Mutasa5_clusts_shorter2_home$perc_time_home)
favstats(Nchelenge5_clusts_shorter2_home$perc_time_home)

perc_time_home2 <- as.data.frame(cbind(c(Choma5_clusts_shorter2_home$perc_time_home, Mutasa5_clusts_shorter2_home$perc_time_home, Nchelenge5_clusts_shorter2_home$perc_time_home),
                                       c(Choma5_clusts_shorter2_home$male, Mutasa5_clusts_shorter2_home$male, Nchelenge5_clusts_shorter2_home$male),
                                       c(rep("Choma",length(Choma5_clusts_shorter2_home$perc_time_home)),
                                         rep("Mutasa",length(Mutasa5_clusts_shorter2_home$perc_time_home)),
                                         rep("Nchelenge",length(Nchelenge5_clusts_shorter2_home$perc_time_home)))))
colnames(perc_time_home2) <- c("values","male","type")
perc_time_home2$values <- as.numeric(as.character(perc_time_home2$values))
perc_time_home2$male <- factor(perc_time_home2$male)
perc_time_home2$type <- factor(perc_time_home2$type)
cols <- brewer.pal(n=4, "Set1")

ggplot(data=perc_time_home2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.3) +
  scale_fill_manual(values = cols2) +
  scale_color_manual(values = cols2) +
  coord_cartesian(xlim=c(-5,110)) +
  scale_x_continuous(breaks=c(seq(0,100,20)))


shapiro.test(perc_time_home2$values) ## not normal
kruskal.test(perc_time_home2$values, perc_time_home2$type) # p = 0.001
pairwise.wilcox.test(perc_time_home2$values, perc_time_home2$type, p.adjust.method = "holm") ## 
#           Choma   Mutasa
# Mutasa    0.008     -     
# Nchelenge  0.003   0.086

my_comp <- list(c("Choma", "Mutasa"), c("Choma","Nchelenge"), c("Mutasa","Nchelenge"))
ggplot(data=perc_time_home2, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = cols2) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.5) +
  stat_compare_means(label.y = 135) 


###### is Choma dry significantly differnet from Nchelenge and Mutasa? ###


perc_time_home2 <- as.data.frame(cbind(c(Choma5_rain_clusts_shorter2_home_dry$perc_time_home, Mutasa5_clusts_shorter2_home$perc_time_home, Nchelenge5_clusts_shorter2_home$perc_time_home),
                                       c(rep("Choma",length(Choma5_rain_clusts_shorter2_home_dry$perc_time_home)),
                                         rep("Mutasa",length(Mutasa5_clusts_shorter2_home$perc_time_home)),
                                         rep("Nchelenge",length(Nchelenge5_clusts_shorter2_home$perc_time_home)))))
colnames(perc_time_home2) <- c("values","type")
perc_time_home2$values <- as.numeric(as.character(perc_time_home2$values))
perc_time_home2$type <- factor(perc_time_home2$type)

shapiro.test(perc_time_home2$values) ## not normal
kruskal.test(perc_time_home2$values, perc_time_home2$type) # p = 0.159
pairwise.wilcox.test(perc_time_home2$values, perc_time_home2$type, p.adjust.method = "holm") ## 
#           Choma   Mutasa
# Mutasa    0.50     -     
# Nchelenge  0.95   0.26


##
##### ??? ####
kruskal.test(perc_time_home2$values[which(perc_time_home2$type=="Choma")], 
             perc_time_home2$male[which(perc_time_home2$type=="Choma")]) ## 
kruskal.test(perc_time_home2$values[which(perc_time_home2$type=="Mutasa")], 
             perc_time_home2$male[which(perc_time_home2$type=="Mutasa")]) ## 
kruskal.test(perc_time_home2$values[which(perc_time_home2$type=="Nchelenge")], 
             perc_time_home2$male[which(perc_time_home2$type=="Nchelenge")]) ## 


ggplot(data=perc_time_home2[which(!(is.na(perc_time_home2$male))),], aes(x=type, y=values, fill=male)) +
  geom_boxplot() +  
  theme_classic() +
  scale_fill_manual("Gender", values = c("red","blue"), labels=c("Female","Male")) +
  scale_x_discrete("Location") + 
  scale_y_continuous("Percent Time at Home") + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))




perc_time_home2$values2 <- perc_time_home2$values/100

descdist(perc_time_home2$values2, boot=1000)
a1 <- fitdist(perc_time_home2$values2, "exp")
a2 <- fitdist(perc_time_home2$values2, "gamma")
a3 <- fitdist(perc_time_home2$values2, "weibull")
a4 <- fitdist(perc_time_home2$values2, "lnorm")
a5 <- fitdist(perc_time_home2$values2, "norm")

denscomp(list(a1,a2,a3,a4,a5)) 
ppcomp(list(a1,a2,a3,a4,a5)) ## none good 



######
#
#### direction ####
temp_Choma <- Choma5_clusts_shorter2_no_home[!duplicated(Choma5_clusts_shorter2_no_home$partid),]
temp_Mutasa <- Mutasa5_clusts_shorter2_no_home[!duplicated(Mutasa5_clusts_shorter2_no_home$partid),]
temp_Mutasa2 <- temp_Mutasa[which(temp_Mutasa$avg_dir_from_home_mag < 75),] ## 164 of 176 parts
temp_Nchelenge <- Nchelenge5_clusts_shorter2_no_home[!duplicated(Nchelenge5_clusts_shorter2_no_home$partid),]
temp_Nchelenge2 <- temp_Nchelenge[which(temp_Nchelenge$avg_dir_from_home_mag < 70),] ## 69 of 71 parts

temp_direction <- as.data.frame(cbind(c(temp_Choma$avg_dir_from_home_bearing,
                                        temp_Mutasa2$avg_dir_from_home_bearing,
                                        temp_Nchelenge2$avg_dir_from_home_bearing),
                                      c(temp_Choma$avg_dir_from_home_mag,
                                        temp_Mutasa2$avg_dir_from_home_mag,
                                        temp_Nchelenge2$avg_dir_from_home_mag),
                                   c(rep("Choma",length(temp_Choma$avg_dir_from_home_bearing)),
                                     rep("Mutasa",length(temp_Mutasa2$avg_dir_from_home_bearing)),
                                     rep("Nchelenge",length(temp_Nchelenge2$avg_dir_from_home_bearing)),
                                     rep("Choma",length(temp_Choma$avg_dir_from_home_bearing)),
                                     rep("Mutasa",length(temp_Mutasa2$avg_dir_from_home_bearing)),
                                     rep("Nchelenge",length(temp_Nchelenge2$avg_dir_from_home_bearing)))))
colnames(temp_direction) <- c("bearing","magnitude","type")
temp_direction$bearing <- as.numeric(as.character(temp_direction$bearing))
temp_direction$magnitude <- as.numeric(as.character(temp_direction$magnitude))
temp_direction$type <- factor(temp_direction$type)
cols <- brewer.pal(n=4, "Set1")

ggplot(data=temp_direction, aes(x=bearing, y=magnitude)) +
  geom_col(aes(fill=type), position=position_dodge2(preserve = "single"), width=3) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual(values=c("magenta","blue","green"))


ggplot(data=temp_direction, aes(x=bearing, y=..count.., fill=type)) +
  geom_histogram(position = "dodge", binwidth = 22.5, alpha=0.8) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~type, nrow=2)+
  scale_fill_manual(values=cols2, name=NULL)



temp_direction <- as.data.frame(cbind(c(temp_Choma$avg_dir_from_home_values,
                                        temp_Mutasa$avg_dir_from_home_values,
                                        temp_Nchelenge$avg_dir_from_home_values),
                                      c(rep("Choma",length(temp_Choma$avg_dir_from_home_values)),
                                        rep("Mutasa",length(temp_Mutasa$avg_dir_from_home_values)),
                                        rep("Nchelenge",length(temp_Nchelenge$avg_dir_from_home_values)))))
colnames(temp_direction) <- c("values","type")
temp_direction$values <- as.numeric(as.character(temp_direction$values))
temp_direction$type <- factor(temp_direction$type)

temp_direction$values[which(temp_direction$values > 348.75)] <- -(360-temp_direction$values[which(temp_direction$values > 348.75)])
ggplot(data=temp_direction, aes(x=values, y=..count.., fill=type)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~type, nrow=2)+
  scale_fill_manual(values=cols2)


temp_direction$values[which(temp_direction$values < 0)] <- 360+temp_direction$values[which(temp_direction$values < 0)]

shapiro.test(temp_direction$values) ## not normal
kruskal.test(temp_direction$values, temp_direction$type) # p = 0.001
pairwise.wilcox.test(temp_direction$values, temp_direction$type, p.adjust.method = "holm") ## 
#           Choma   Mutasa
# Mutasa    0.002     -     
# Nchelenge  0.002   0.576

my_comp <- list(c("Choma", "Mutasa"), c("Choma","Nchelenge"), c("Mutasa","Nchelenge"))
ggplot(data=temp_direction, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = cols2) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.5) +
  stat_compare_means(label.y = 480) 





## average direction calc
temp_direction$values2 <- 90 - temp_direction$values
temp_direction$values2[which(temp_direction$values2 < 0)] <- temp_direction$values2[which(temp_direction$values2 < 0)] + 360
temp_direction$values3 <- deg2rad(temp_direction$values2)

temp_direction$avg_value <- as.numeric(NA)
for(ii in 1:nlevels(temp_direction$type)){
  level <- levels(temp_direction$type)[ii]
  temp <- temp_direction[which(temp_direction$type == level),]
  sum_cos_dist <- sum(cos(temp$values3))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$values3))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  temp_direction$avg_value[which(temp_direction$type == level)] <- avg_deg
}

temp_direction$avg_value2 <- 90 - temp_direction$avg_value
temp_direction$avg_value2[which(temp_direction$avg_value2 < 0)] <- temp_direction$avg_value2[which(temp_direction$avg_value2 < 0)] + 360
favstats(temp_direction$avg_value2~temp_direction$type)
## Choma: 62.75 --> ENE
## Mutasa: 276.47 --> W
## Nchelenge: 274.12 --> W


#
#
#
### Time at locations ####
favstats(Choma5_clusts_shorter2_no_home$clust_trips_time_new*24)
clusts_median_Choma <- favstats((Choma5_clusts_shorter2_no_home$clust_trips_time_new*24)~Choma5_clusts_shorter2_no_home$partid)$median
favstats(clusts_median_Choma)
favstats(Mutasa5_clusts_shorter2_no_home$clust_trips_time_new*24)
clusts_median_Mutasa <- favstats((Mutasa5_clusts_shorter2_no_home$clust_trips_time_new*24)~Mutasa5_clusts_shorter2_no_home$partid)$median
favstats(clusts_median_Mutasa)
clusts_median_Mutasa2 <- clusts_median_Mutasa[which(clusts_median_Mutasa < 150)]
favstats(clusts_median_Mutasa2)
favstats(Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new*24)
clusts_median_Nchelenge <- favstats((Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new*24)~Nchelenge5_clusts_shorter2_no_home$partid)$median
favstats(clusts_median_Nchelenge)
#
#

clust_time2 <- as.data.frame(cbind(c(clusts_median_Choma, clusts_median_Mutasa2, clusts_median_Nchelenge),
                                   c(rep("Choma",length(clusts_median_Choma)),
                                     rep("Mutasa",length(clusts_median_Mutasa2)),
                                     rep("Nchelenge",length(clusts_median_Nchelenge)))))
colnames(clust_time2) <- c("values","type")
clust_time2$values <- as.numeric(as.character(clust_time2$values))
clust_time2$type <- factor(clust_time2$type)
cols <- brewer.pal(n=4, "Set1")

ggplot(data=clust_time2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=12) + 
  scale_fill_manual(values = cols2) +
  scale_color_manual(values = cols2) +
  coord_cartesian(xlim=c(-3,113)) +
  scale_x_continuous(breaks=c(seq(0,200,10)))



shapiro.test(clust_time2$values) ## not normal
kruskal.test(clust_time2$values, clust_time2$type) # p < 2e-16
pairwise.wilcox.test(clust_time2$values, clust_time2$type, p.adjust.method = "holm") ## 
#           Choma   Mutasa
# Mutasa    < 2e-16     -     
# Nchelenge  0.02    7e-15

my_comp <- list(c("Choma", "Mutasa"), c("Choma","Nchelenge"), c("Mutasa","Nchelenge"))
ggplot(data=clust_time2, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = cols2) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.5) +
  stat_compare_means(label.y = 155) 
###
#### Biting time metrics #####
### number locations ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home_biting_places2$clust~Choma5_clusts_shorter2_no_home_biting_places2$partid)$n
only_home_Choma <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2$partid))
loc_counts_Choma2 <- c(loc_counts_Choma,rep(0, only_home_Choma))

loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home_biting_places2$clust~Mutasa5_clusts_shorter2_no_home_biting_places2$partid)$n
only_home_Mutasa <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2$partid))
loc_counts_Mutasa2 <- c(loc_counts_Mutasa,rep(0, only_home_Mutasa))

loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2$partid)$n
only_home_Nchelenge <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2$partid))
loc_counts_Nchelenge2 <- c(loc_counts_Nchelenge,rep(0, only_home_Nchelenge))

locs_count2 <- as.data.frame(cbind(c(loc_counts_Choma2, loc_counts_Mutasa2, loc_counts_Nchelenge2),
                                   c(rep("Choma",length(loc_counts_Choma2)),
                                     rep("Mutasa",length(loc_counts_Mutasa2)),
                                     rep("Nchelenge",length(loc_counts_Nchelenge2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)
cols <- brewer.pal(n=4, "Set1")

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = cols2) +
  scale_color_manual(values = cols2) +
  coord_cartesian(xlim=c(-0.75,13)) +
  scale_x_continuous(breaks=c(seq(0,14,2)))
##
### time per location ####
clusts_median_Choma <- (favstats(Choma5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2$partid)$median)*24
clusts_median_Mutasa <- (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2$partid)$median)*24
clusts_median_Nchelenge <- (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2$partid)$median)*24
clusts_median_Mutasa2 <- clusts_median_Mutasa[which(clusts_median_Mutasa < 500)]

clust_time2 <- as.data.frame(cbind(c(clusts_median_Choma, clusts_median_Mutasa2, clusts_median_Nchelenge),
                                   c(rep("Choma",length(clusts_median_Choma)),
                                     rep("Mutasa",length(clusts_median_Mutasa2)),
                                     rep("Nchelenge",length(clusts_median_Nchelenge)))))
colnames(clust_time2) <- c("values","type")
clust_time2$values <- as.numeric(as.character(clust_time2$values))
clust_time2$type <- factor(clust_time2$type)
cols <- brewer.pal(n=4, "Set1")

ggplot(data=clust_time2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 50, alpha=0.8) +
  geom_density(aes(color=type), adjust=5) + 
  scale_fill_manual(values = cols2) +
  scale_color_manual(values = cols2) +
  coord_cartesian(xlim=c(-18,300)) +
  scale_x_continuous(breaks=c(seq(0,300,50)))
#
### number trips ####
### PER PART ###
trip_count_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3$partid)$median
trip_count_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3$partid)$median
trip_count_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3$partid)$median
trip_count_median_Choma2 <- trip_count_median_Choma[which(trip_count_median_Choma < 20)]

trips_all <- as.data.frame(cbind(c(trip_count_median_Choma2, trip_count_median_Mutasa, trip_count_median_Nchelenge),
                                 c(rep("Choma",length(trip_count_median_Choma2)),
                                   rep("Mutasa",length(trip_count_median_Mutasa)),
                                   rep("Nchelenge",length(trip_count_median_Nchelenge)))))
colnames(trips_all) <- c("values","type")
trips_all$values <- as.numeric(as.character(trips_all$values))
trips_all$type <- factor(trips_all$type)
cols <- brewer.pal(n=4, "Set1")

ggplot(data=trips_all, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=4) +
  scale_fill_manual(values = cols2) +
  scale_color_manual(values = cols2) +
  coord_cartesian(xlim=c(0,13)) +
  scale_x_continuous(breaks=c(seq(0,30,2)))
#

trips_all2 <- as.data.frame(cbind(c(Choma5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2, 
                                   Mutasa5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2, 
                                   Nchelenge5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2),
                                 c(rep("Choma",length(Choma5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2)),
                                   rep("Mutasa",length(Mutasa5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2)),
                                   rep("Nchelenge",length(Nchelenge5_clust_trips_shorter2_no_home_biting_places3$trip_count_new2)))))
colnames(trips_all2) <- c("values","type")
trips_all2$values <- as.numeric(as.character(trips_all2$values))
trips_all2$type <- factor(trips_all2$type)

table(trips_all2$type, as.factor(trips_all2$values))

#
##### percent biting time spent at home vs elsewhere #####
perc_time_home2 <- as.data.frame(cbind(c(Choma5_clusts_shorter2_biting_places3b$percent_home_biting_time, 
                                         Mutasa5_clusts_shorter2_biting_places3b$percent_home_biting_time, 
                                         Nchelenge5_clusts_shorter2_biting_places3b$percent_home_biting_time),
                                       c(rep("Choma",length(Choma5_clusts_shorter2_biting_places3b$percent_home_biting_time)),
                                         rep("Mutasa",length(Mutasa5_clusts_shorter2_biting_places3b$percent_home_biting_time)),
                                         rep("Nchelenge",length(Nchelenge5_clusts_shorter2_biting_places3b$percent_home_biting_time)))))
colnames(perc_time_home2) <- c("values","type")
perc_time_home2$values <- as.numeric(as.character(perc_time_home2$values))
perc_time_home2$type <- factor(perc_time_home2$type)
cols <- brewer.pal(n=4, "Set1")

ggplot(data=perc_time_home2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = cols2) +
  scale_color_manual(values = cols2) +
  coord_cartesian(xlim=c(-5,110)) +
  scale_x_continuous(breaks=c(seq(0,100,20)))
##
###############
### number of locations vs number of trips ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,9)]

trips_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home3$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]
trips_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home3$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]
trips_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home3$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]


locs_and_trips_count <- as.data.frame(cbind(c(loc_counts_Choma$`Choma5_clusts_shorter2_no_home$partid`, 
                                     loc_counts_Mutasa$`Mutasa5_clusts_shorter2_no_home$partid`, 
                                     loc_counts_Nchelenge$`Nchelenge5_clusts_shorter2_no_home$partid`),
                                   c(loc_counts_Choma$n, 
                                     loc_counts_Mutasa$n, 
                                     loc_counts_Nchelenge$n),
                                   c(trips_part_median_Choma$median, 
                                     trips_part_median_Mutasa$median, 
                                     trips_part_median_Nchelenge$median),
                                   c(rep("Choma",length(loc_counts_Choma$n)),
                                     rep("Mutasa",length(loc_counts_Mutasa$n)),
                                     rep("Nchelenge",length(loc_counts_Nchelenge$n)))))
colnames(locs_and_trips_count) <- c("partid","num_locs","num_trips","loc")
locs_and_trips_count$num_locs <- as.numeric(as.character(locs_and_trips_count$num_locs))
locs_and_trips_count$num_trips <- as.numeric(as.character(locs_and_trips_count$num_trips))
locs_and_trips_count$partid <- factor(locs_and_trips_count$partid)
locs_and_trips_count$loc <- factor(locs_and_trips_count$loc)

temp <- as.data.frame(cbind(c(Choma5_clusts_shorter2_home$partid, 
                              Mutasa5_clusts_shorter2_home$partid, 
                              Nchelenge5_clusts_shorter2_home$partid),
                            c(Choma5_clusts_shorter2_home$perc_time_home, 
                              Mutasa5_clusts_shorter2_home$perc_time_home, 
                              Nchelenge5_clusts_shorter2_home$perc_time_home)))
colnames(temp) <- c("partid","perc_time_home")
temp$perc_time_home <- as.numeric(as.character(temp$perc_time_home))
temp$partid <- factor(temp$partid)

locs_and_trips_count2 <- droplevels(merge(locs_and_trips_count, temp, by="partid"))
locs_and_trips_count2$time_out <- (100-locs_and_trips_count2$perc_time_home)/100

ggplot(data=locs_and_trips_count2, aes(x=num_locs, y=num_trips, fill=loc, size=time_out)) +
  geom_point(position=position_jitter(width=0.25, height=0.25), shape=21, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,15,1)))


ggplot(data=locs_and_trips_count2, aes(x=num_locs, y=time_out, fill=loc)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,1,0.1)))

ggplot(data=locs_and_trips_count2, aes(x=num_trips, y=time_out, fill=loc)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,15,1))) +
  scale_y_continuous(breaks=c(seq(0,1,0.1)))

#



### number of locations vs time per trip ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,9)]

trips_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home3$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]
trips_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home3$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]
trips_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home3$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]

trip_time_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home$partid)[,c(1,4)]
trip_time_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home$partid)[,c(1,4)]
trip_time_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home$partid)[,c(1,4)]

locs_and_trips_time <- as.data.frame(cbind(c(loc_counts_Choma$`Choma5_clusts_shorter2_no_home$partid`, 
                                              loc_counts_Mutasa$`Mutasa5_clusts_shorter2_no_home$partid`, 
                                              loc_counts_Nchelenge$`Nchelenge5_clusts_shorter2_no_home$partid`),
                                            c(loc_counts_Choma$n, 
                                              loc_counts_Mutasa$n, 
                                              loc_counts_Nchelenge$n),
                                            c(trip_time_part_median_Choma$median, 
                                              trip_time_part_median_Mutasa$median, 
                                              trip_time_part_median_Nchelenge$median),
                                            c(rep("Choma",length(loc_counts_Choma$n)),
                                              rep("Mutasa",length(loc_counts_Mutasa$n)),
                                              rep("Nchelenge",length(loc_counts_Nchelenge$n)))))
colnames(locs_and_trips_time) <- c("partid","num_locs","trip_time","loc")
locs_and_trips_time$num_locs <- as.numeric(as.character(locs_and_trips_time$num_locs))
locs_and_trips_time$trip_time <- as.numeric(as.character(locs_and_trips_time$trip_time))
locs_and_trips_time$partid <- factor(locs_and_trips_time$partid)
locs_and_trips_time$loc <- factor(locs_and_trips_time$loc)
## remove outliers ##
locs_and_trips_time2 <- locs_and_trips_time[which(locs_and_trips_time$trip_time < 30),]


temp <- as.data.frame(cbind(c(Choma5_clusts_shorter2_home$partid, 
                              Mutasa5_clusts_shorter2_home$partid, 
                              Nchelenge5_clusts_shorter2_home$partid),
                            c(Choma5_clusts_shorter2_home$perc_time_home, 
                              Mutasa5_clusts_shorter2_home$perc_time_home, 
                              Nchelenge5_clusts_shorter2_home$perc_time_home)))
colnames(temp) <- c("partid","perc_time_home")
temp$perc_time_home <- as.numeric(as.character(temp$perc_time_home))
temp$partid <- factor(temp$partid)

locs_and_trips_time2b <- droplevels(merge(locs_and_trips_time2, temp, by="partid"))
locs_and_trips_time2b$time_out <- (100-locs_and_trips_time2b$perc_time_home)/100

cvdPlot(ggplot(data=locs_and_trips_time2b, aes(x=num_locs, y=trip_time, fill=loc, size=time_out)) +
  geom_point(position=position_jitter(width=0.25, height=0), shape=21, col="black") +
  scale_fill_manual("Site", values = c(SteppedSequential5Steps[c(21,18,15)])) +
  scale_size("Percent time\noutside home", breaks=c(0.25,0.50,0.75),labels=c("25%","50%","75%")) +
  scale_x_continuous("Number of locations visited", breaks=c(seq(0,26,2))) +
  scale_y_continuous("Meidan time per trip", breaks=c(seq(0,13,1))) +
  theme_classic() +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=12),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12)) +
  guides(fill = guide_legend(override.aes = list(size=5))), layout = "desaturate")


ggplot(data=locs_and_trips_time2b, aes(x=num_locs, y=trip_time, fill=loc, size=time_out)) +
  geom_point(position=position_jitter(width=0.2, height=0), shape=21, col="black") +
  scale_fill_manual("Site", values = c(SteppedSequential5Steps[c(21,18,15)])) +
  scale_size("Percent time\noutside home", breaks=c(0.25,0.50,0.75),labels=c("25%","50%","75%")) +
  scale_x_continuous("Number of locations visited", breaks=c(seq(0,26,2))) +
  scale_y_continuous("Meidan time per trip (hr)", breaks=c(seq(0,13,1))) +
  theme_classic() +
  theme(axis.title=element_text(size=16),
        axis.text = element_text(size=14),
        legend.title=element_text(size=16),
        legend.text=element_text(size=14)) +
  guides(fill = guide_legend(override.aes = list(size=5)))


ggplot(data=locs_and_trips_time2b, aes(x=trip_time, y=time_out, fill=loc)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,13,1))) +
  scale_y_continuous(breaks=c(seq(0,1,0.1)))
#


locs_and_trips_count_and_trips_time <- merge(locs_and_trips_count2, locs_and_trips_time2b[,c(1,3)], by="partid")

ggplot(data=locs_and_trips_count_and_trips_time, aes(x=num_trips, y=trip_time, fill=loc)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_fill_manual(values = cols2) +
  scale_y_continuous(breaks=c(seq(0,13,1))) +
  scale_x_continuous(breaks=c(seq(0,15,1)))
#


# plot_ly(x=locs_and_trips_count_and_trips_time$num_locs, 
#         y=locs_and_trips_count_and_trips_time$num_trips, 
#         z=locs_and_trips_count_and_trips_time$trip_time, 
#         type="scatter3d", 
#         color=locs_and_trips_count_and_trips_time$loc,
#         size=locs_and_trips_count_and_trips_time$time_out, sizes=c(10,400))
# 
# 
# plot_ly(x=locs_and_trips_count_and_trips_time$num_locs, 
#         y=locs_and_trips_count_and_trips_time$trip_time, 
#         z=locs_and_trips_count_and_trips_time$time_out,
#         type="scatter3d", 
#         color=locs_and_trips_count_and_trips_time$loc)
# ##

### number of locations vs distance ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,9)]

loc_dist_median_Choma <- favstats(Choma5_clusts_shorter2_no_home$hhdist_km_haversine ~ Choma5_clusts_shorter2_no_home$partid)[,c(1,4)]
loc_dist_median_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$hhdist_km_haversine ~ Mutasa5_clusts_shorter2_no_home$partid)[,c(1,4)]
loc_dist_median_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine ~ Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,4)]

locs_and_dist <- as.data.frame(cbind(c(loc_counts_Choma$`Choma5_clusts_shorter2_no_home$partid`, 
                                             loc_counts_Mutasa$`Mutasa5_clusts_shorter2_no_home$partid`, 
                                             loc_counts_Nchelenge$`Nchelenge5_clusts_shorter2_no_home$partid`),
                                           c(loc_counts_Choma$n, 
                                             loc_counts_Mutasa$n, 
                                             loc_counts_Nchelenge$n),
                                           c(loc_dist_median_Choma$median, 
                                             loc_dist_median_Mutasa$median, 
                                             loc_dist_median_Nchelenge$median),
                                           c(rep("Choma",length(loc_counts_Choma$n)),
                                             rep("Mutasa",length(loc_counts_Mutasa$n)),
                                             rep("Nchelenge",length(loc_counts_Nchelenge$n)))))
colnames(locs_and_dist) <- c("partid","num_locs","loc_dist","loc")
locs_and_dist$num_locs <- as.numeric(as.character(locs_and_dist$num_locs))
locs_and_dist$loc_dist <- as.numeric(as.character(locs_and_dist$loc_dist))
locs_and_dist$partid <- factor(locs_and_dist$partid)
locs_and_dist$loc <- factor(locs_and_dist$loc)
## remove outliers ##
locs_and_dist2 <- locs_and_dist[which(locs_and_dist$loc_dist < 85),]

ggplot(data=locs_and_dist2, aes(x=num_locs, y=loc_dist, fill=loc)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))

#
##
#### time per trips vs distance ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,9)]

loc_dist_median_Choma <- favstats(Choma5_clusts_shorter2_no_home$hhdist_km_haversine ~ Choma5_clusts_shorter2_no_home$partid)[,c(1,4)]
loc_dist_median_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$hhdist_km_haversine ~ Mutasa5_clusts_shorter2_no_home$partid)[,c(1,4)]
loc_dist_median_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine ~ Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,4)]

trip_time_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home$partid)[,c(1,4)]
trip_time_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home$partid)[,c(1,4)]
trip_time_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home$partid)[,c(1,4)]

locs_num_dist_and_trip_time <- as.data.frame(cbind(c(loc_counts_Choma$`Choma5_clusts_shorter2_no_home$partid`, 
                                       loc_counts_Mutasa$`Mutasa5_clusts_shorter2_no_home$partid`, 
                                       loc_counts_Nchelenge$`Nchelenge5_clusts_shorter2_no_home$partid`),
                                     c(loc_counts_Choma$n, 
                                       loc_counts_Mutasa$n, 
                                       loc_counts_Nchelenge$n),
                                     c(loc_dist_median_Choma$median, 
                                       loc_dist_median_Mutasa$median, 
                                       loc_dist_median_Nchelenge$median),
                                     c(trip_time_part_median_Choma$median, 
                                       trip_time_part_median_Mutasa$median, 
                                       trip_time_part_median_Nchelenge$median),
                                     c(rep("Choma",length(loc_counts_Choma$n)),
                                       rep("Mutasa",length(loc_counts_Mutasa$n)),
                                       rep("Nchelenge",length(loc_counts_Nchelenge$n)))))
colnames(locs_num_dist_and_trip_time) <- c("partid","num_locs","loc_dist","trip_time", "loc")
locs_num_dist_and_trip_time$num_locs <- as.numeric(as.character(locs_num_dist_and_trip_time$num_locs))
locs_num_dist_and_trip_time$loc_dist <- as.numeric(as.character(locs_num_dist_and_trip_time$loc_dist))
locs_num_dist_and_trip_time$trip_time <- as.numeric(as.character(locs_num_dist_and_trip_time$trip_time))
locs_num_dist_and_trip_time$partid <- factor(locs_num_dist_and_trip_time$partid)
locs_num_dist_and_trip_time$loc <- factor(locs_num_dist_and_trip_time$loc)
## remove outliers ##
locs_num_dist_and_trip_time2 <- locs_num_dist_and_trip_time[which(locs_num_dist_and_trip_time$loc_dist < 85),]
locs_num_dist_and_trip_time2b <- locs_num_dist_and_trip_time2[which(locs_num_dist_and_trip_time2$trip_time < 30),]


ggplot(data=locs_num_dist_and_trip_time2b, aes(x=num_locs, y=loc_dist, fill=loc)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))


ggplot(data=locs_num_dist_and_trip_time2b, aes(x=trip_time, y=loc_dist, fill=loc)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,13,1))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))


###
###############
### number of locations vs average number of trips ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,9)]

trips_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home3$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]
trips_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home3$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]
trips_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home3$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]


locs_and_trips_count <- as.data.frame(cbind(c(loc_counts_Choma$`Choma5_clusts_shorter2_no_home$partid`, 
                                              loc_counts_Mutasa$`Mutasa5_clusts_shorter2_no_home$partid`, 
                                              loc_counts_Nchelenge$`Nchelenge5_clusts_shorter2_no_home$partid`),
                                            c(loc_counts_Choma$n, 
                                              loc_counts_Mutasa$n, 
                                              loc_counts_Nchelenge$n),
                                            c(trips_part_median_Choma$median, 
                                              trips_part_median_Mutasa$median, 
                                              trips_part_median_Nchelenge$median),
                                            c(rep("Choma",length(loc_counts_Choma$n)),
                                              rep("Mutasa",length(loc_counts_Mutasa$n)),
                                              rep("Nchelenge",length(loc_counts_Nchelenge$n)))))
colnames(locs_and_trips_count) <- c("partid","num_locs","num_trips","loc")
locs_and_trips_count$num_locs <- as.numeric(as.character(locs_and_trips_count$num_locs))
locs_and_trips_count$num_trips <- as.numeric(as.character(locs_and_trips_count$num_trips))
locs_and_trips_count$partid <- factor(locs_and_trips_count$partid)
locs_and_trips_count$loc <- factor(locs_and_trips_count$loc)

temp <- as.data.frame(cbind(c(Choma5_clusts_shorter2_home$partid, 
                              Mutasa5_clusts_shorter2_home$partid, 
                              Nchelenge5_clusts_shorter2_home$partid),
                            c(Choma5_clusts_shorter2_home$perc_time_home, 
                              Mutasa5_clusts_shorter2_home$perc_time_home, 
                              Nchelenge5_clusts_shorter2_home$perc_time_home)))
colnames(temp) <- c("partid","perc_time_home")
temp$perc_time_home <- as.numeric(as.character(temp$perc_time_home))
temp$partid <- factor(temp$partid)

locs_and_trips_count2 <- droplevels(merge(locs_and_trips_count, temp, by="partid"))
locs_and_trips_count2$time_out <- (100-locs_and_trips_count2$perc_time_home)/100

ggplot(data=locs_and_trips_count2, aes(x=num_locs, y=num_trips, fill=loc, size=time_out)) +
  geom_point(position=position_jitter(width=0.25, height=0.25), shape=21, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,15,1)))

ggplot(data=locs_and_trips_count2[which(locs_and_trips_count2$loc == "Choma"),], aes(x=num_locs, y=num_trips, size=time_out)) +
  geom_point(position=position_jitter(width=0.25, height=0.25), shape=21, col="black") +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,15,1)))

ggplot(data=locs_and_trips_count2[which(locs_and_trips_count2$loc == "Mutasa"),], aes(x=num_locs, y=num_trips, size=time_out)) +
  geom_point(position=position_jitter(width=0.25, height=0.25), shape=21, col="black") +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,15,1)))

ggplot(data=locs_and_trips_count2[which(locs_and_trips_count2$loc == "Nchelenge"),], aes(x=num_locs, y=num_trips, size=time_out)) +
  geom_point(position=position_jitter(width=0.25, height=0.25), shape=21, col="black") +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,15,1)))
###
### number of locations vs time per trip ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,9)]

trips_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home3$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]
trips_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home3$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]
trips_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home3$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]

trip_time_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home$partid)[,c(1,4)]
trip_time_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home$partid)[,c(1,4)]
trip_time_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home$partid)[,c(1,4)]

locs_and_trips_time <- as.data.frame(cbind(c(loc_counts_Choma$`Choma5_clusts_shorter2_no_home$partid`, 
                                             loc_counts_Mutasa$`Mutasa5_clusts_shorter2_no_home$partid`, 
                                             loc_counts_Nchelenge$`Nchelenge5_clusts_shorter2_no_home$partid`),
                                           c(loc_counts_Choma$n, 
                                             loc_counts_Mutasa$n, 
                                             loc_counts_Nchelenge$n),
                                           c(trip_time_part_median_Choma$median, 
                                             trip_time_part_median_Mutasa$median, 
                                             trip_time_part_median_Nchelenge$median),
                                           c(rep("Choma",length(loc_counts_Choma$n)),
                                             rep("Mutasa",length(loc_counts_Mutasa$n)),
                                             rep("Nchelenge",length(loc_counts_Nchelenge$n)))))
colnames(locs_and_trips_time) <- c("partid","num_locs","trip_time","loc")
locs_and_trips_time$num_locs <- as.numeric(as.character(locs_and_trips_time$num_locs))
locs_and_trips_time$trip_time <- as.numeric(as.character(locs_and_trips_time$trip_time))
locs_and_trips_time$partid <- factor(locs_and_trips_time$partid)
locs_and_trips_time$loc <- factor(locs_and_trips_time$loc)
## remove outliers ##
locs_and_trips_time2 <- locs_and_trips_time[which(locs_and_trips_time$trip_time < 30),]


temp <- as.data.frame(cbind(c(Choma5_clusts_shorter2_home$partid, 
                              Mutasa5_clusts_shorter2_home$partid, 
                              Nchelenge5_clusts_shorter2_home$partid),
                            c(Choma5_clusts_shorter2_home$perc_time_home, 
                              Mutasa5_clusts_shorter2_home$perc_time_home, 
                              Nchelenge5_clusts_shorter2_home$perc_time_home)))
colnames(temp) <- c("partid","perc_time_home")
temp$perc_time_home <- as.numeric(as.character(temp$perc_time_home))
temp$partid <- factor(temp$partid)

locs_and_trips_time2b <- droplevels(merge(locs_and_trips_time2, temp, by="partid"))
locs_and_trips_time2b$time_out <- (100-locs_and_trips_time2b$perc_time_home)/100
locs_and_trips_time2b$time_in <- (locs_and_trips_time2b$perc_time_home)/100

ggplot(data=locs_and_trips_time2b, aes(x=num_locs, y=trip_time, fill=loc, size=time_out)) +
  geom_point(position=position_jitter(width=0.25, height=0), shape=21, col="black") +
  scale_fill_manual(values = cols2) +
  scale_size_area(max_size = 7) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,13,1)))


ggplot(data=locs_and_trips_time2b[which(locs_and_trips_time2b$loc == "Choma"),], aes(x=num_locs, y=trip_time, size=time_out)) +
  geom_point(position=position_jitter(width=0.25, height=0), shape=21, col="black",fill=cols[4]) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,13,1)))

ggplot(data=locs_and_trips_time2b[which(locs_and_trips_time2b$loc == "Mutasa"),], aes(x=num_locs, y=trip_time, size=time_out)) +
  geom_point(position=position_jitter(width=0.25, height=0), shape=21, col="black",fill=cols[2]) +
  scale_size_area(max_size = 6) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,13,1)))

ggplot(data=locs_and_trips_time2b[which(locs_and_trips_time2b$loc == "Nchelenge"),], aes(x=num_locs, y=trip_time, size=time_out)) +
  geom_point(position=position_jitter(width=0.25, height=0), shape=21, col="black",fill=cols[3]) +
  scale_size_area(max_size = 6) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,13,1)))


ggplot(data=locs_and_trips_time2b, aes(x=num_locs, y=trip_time, fill=loc)) +
  geom_point(position=position_jitter(width=0.25, height=0), shape=21, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,13,1)))

ggplot(data=locs_and_trips_time2b, aes(x=num_locs, y=trip_time)) +
  geom_point(position=position_jitter(width=0.25, height=0), shape=21, col="black") +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,13,1)))


aa <- locs_and_trips_time2b[,c(2,3,6)]
aa2 <- cor(aa)
corrplot(aa2,method="number",type="lower")
rcorr(as.matrix(aa), type="pearson")
chart.Correlation(aa, histogram=TRUE, pch=19, method="pearson")
##


ggplot(data=locs_and_trips_time2b, aes(x=num_locs, y=trip_time, fill=loc, size=time_out)) +
  geom_point(position=position_jitter(width=0.25, height=0), shape=21, col="black") +
  scale_fill_manual(values = cols2) +
  scale_size_area(max_size = 7) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,13,1)))

ggplot(data=locs_and_trips_time2b, aes(x=num_locs, y=time_out, fill=loc, size=trip_time)) +
  geom_point(position=position_jitter(width=0.25, height=0), shape=21, col="black") +
  scale_fill_manual(values = cols2) +
  scale_size_area(max_size = 7) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,1,0.1)))

bb <- locs_and_trips_time2b[which(locs_and_trips_time2b$loc=="Choma"),c(2,3,6)]
rcorr(as.matrix(bb), type="pearson")
chart.Correlation(bb, histogram=TRUE, pch=19, method="pearson")

cc <- locs_and_trips_time2b[which(locs_and_trips_time2b$loc=="Mutasa"),c(2,3,6)]
rcorr(as.matrix(cc), type="pearson")
chart.Correlation(cc, histogram=TRUE, pch=19, method="pearson")

dd <- locs_and_trips_time2b[which(locs_and_trips_time2b$loc=="Nchelenge"),c(2,3,6)]
rcorr(as.matrix(dd), type="pearson")
chart.Correlation(dd, histogram=TRUE, pch=19, method="pearson")


##

### number of locations vs distance ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,9)]

loc_dist_median_Choma <- favstats(Choma5_clusts_shorter2_no_home$hhdist_km_haversine ~ Choma5_clusts_shorter2_no_home$partid)[,c(1,4)]
loc_dist_median_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$hhdist_km_haversine ~ Mutasa5_clusts_shorter2_no_home$partid)[,c(1,4)]
loc_dist_median_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine ~ Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,4)]

locs_and_dist <- as.data.frame(cbind(c(loc_counts_Choma$`Choma5_clusts_shorter2_no_home$partid`, 
                                       loc_counts_Mutasa$`Mutasa5_clusts_shorter2_no_home$partid`, 
                                       loc_counts_Nchelenge$`Nchelenge5_clusts_shorter2_no_home$partid`),
                                     c(loc_counts_Choma$n, 
                                       loc_counts_Mutasa$n, 
                                       loc_counts_Nchelenge$n),
                                     c(loc_dist_median_Choma$median, 
                                       loc_dist_median_Mutasa$median, 
                                       loc_dist_median_Nchelenge$median),
                                     c(rep("Choma",length(loc_counts_Choma$n)),
                                       rep("Mutasa",length(loc_counts_Mutasa$n)),
                                       rep("Nchelenge",length(loc_counts_Nchelenge$n)))))
colnames(locs_and_dist) <- c("partid","num_locs","loc_dist","loc")
locs_and_dist$num_locs <- as.numeric(as.character(locs_and_dist$num_locs))
locs_and_dist$loc_dist <- as.numeric(as.character(locs_and_dist$loc_dist))
locs_and_dist$partid <- factor(locs_and_dist$partid)
locs_and_dist$loc <- factor(locs_and_dist$loc)
## remove outliers ##
locs_and_dist2 <- locs_and_dist[which(locs_and_dist$loc_dist < 85),]

ggplot(data=locs_and_dist2, aes(x=num_locs, y=loc_dist, fill=loc)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))

ggplot(data=locs_and_dist2[which(locs_and_dist2$loc == "Choma"),], aes(x=num_locs, y=loc_dist)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))

ggplot(data=locs_and_dist2[which(locs_and_dist2$loc == "Mutasa"),], aes(x=num_locs, y=loc_dist)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))

ggplot(data=locs_and_dist2[which(locs_and_dist2$loc == "Nchelenge"),], aes(x=num_locs, y=loc_dist)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))

#
##
#### time per trips vs distance ####
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,9)]

loc_dist_median_Choma <- favstats(Choma5_clusts_shorter2_no_home$hhdist_km_haversine ~ Choma5_clusts_shorter2_no_home$partid)[,c(1,4)]
loc_dist_median_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$hhdist_km_haversine ~ Mutasa5_clusts_shorter2_no_home$partid)[,c(1,4)]
loc_dist_median_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine ~ Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,4)]

trip_time_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home$partid)[,c(1,4)]
trip_time_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home$partid)[,c(1,4)]
trip_time_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home$partid)[,c(1,4)]

locs_num_dist_and_trip_time <- as.data.frame(cbind(c(loc_counts_Choma$`Choma5_clusts_shorter2_no_home$partid`, 
                                                     loc_counts_Mutasa$`Mutasa5_clusts_shorter2_no_home$partid`, 
                                                     loc_counts_Nchelenge$`Nchelenge5_clusts_shorter2_no_home$partid`),
                                                   c(loc_counts_Choma$n, 
                                                     loc_counts_Mutasa$n, 
                                                     loc_counts_Nchelenge$n),
                                                   c(loc_dist_median_Choma$median, 
                                                     loc_dist_median_Mutasa$median, 
                                                     loc_dist_median_Nchelenge$median),
                                                   c(trip_time_part_median_Choma$median, 
                                                     trip_time_part_median_Mutasa$median, 
                                                     trip_time_part_median_Nchelenge$median),
                                                   c(rep("Choma",length(loc_counts_Choma$n)),
                                                     rep("Mutasa",length(loc_counts_Mutasa$n)),
                                                     rep("Nchelenge",length(loc_counts_Nchelenge$n)))))
colnames(locs_num_dist_and_trip_time) <- c("partid","num_locs","loc_dist","trip_time", "loc")
locs_num_dist_and_trip_time$num_locs <- as.numeric(as.character(locs_num_dist_and_trip_time$num_locs))
locs_num_dist_and_trip_time$loc_dist <- as.numeric(as.character(locs_num_dist_and_trip_time$loc_dist))
locs_num_dist_and_trip_time$trip_time <- as.numeric(as.character(locs_num_dist_and_trip_time$trip_time))
locs_num_dist_and_trip_time$partid <- factor(locs_num_dist_and_trip_time$partid)
locs_num_dist_and_trip_time$loc <- factor(locs_num_dist_and_trip_time$loc)
## remove outliers ##
locs_num_dist_and_trip_time2 <- locs_num_dist_and_trip_time[which(locs_num_dist_and_trip_time$loc_dist < 85),]
locs_num_dist_and_trip_time2b <- locs_num_dist_and_trip_time2[which(locs_num_dist_and_trip_time2$trip_time < 30),]


ggplot(data=locs_num_dist_and_trip_time2b, aes(x=num_locs, y=loc_dist, fill=loc)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))

ggplot(data=locs_num_dist_and_trip_time2b[which(locs_num_dist_and_trip_time2b$loc == "Choma"),], aes(x=num_locs, y=loc_dist)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))

ggplot(data=locs_num_dist_and_trip_time2b[which(locs_num_dist_and_trip_time2b$loc == "Mutasa"),], aes(x=num_locs, y=loc_dist)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))

ggplot(data=locs_num_dist_and_trip_time2b[which(locs_num_dist_and_trip_time2b$loc == "Nchelenge"),], aes(x=num_locs, y=loc_dist)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_x_continuous(breaks=c(seq(0,26,2))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))






ggplot(data=locs_num_dist_and_trip_time2b, aes(x=trip_time, y=loc_dist, fill=loc)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_fill_manual(values = cols2) +
  scale_x_continuous(breaks=c(seq(0,13,1))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))

ggplot(data=locs_num_dist_and_trip_time2b[which(locs_num_dist_and_trip_time2b$loc == "Choma"),], aes(x=trip_time, y=loc_dist)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_x_continuous(breaks=c(seq(0,13,1))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))

ggplot(data=locs_num_dist_and_trip_time2b[which(locs_num_dist_and_trip_time2b$loc == "Mutasa"),], aes(x=trip_time, y=loc_dist)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_x_continuous(breaks=c(seq(0,13,1))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))

ggplot(data=locs_num_dist_and_trip_time2b[which(locs_num_dist_and_trip_time2b$loc == "Nchelenge"),], aes(x=trip_time, y=loc_dist)) +
  geom_point(position=position_jitter(width=0.3, height=0), shape=21, size=2, col="black") +
  scale_x_continuous(breaks=c(seq(0,13,1))) +
  scale_y_continuous(breaks=c(seq(0,85,10)))



###

###############
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
############ ############ ############ ############ ############ ############ ############ 
#############################################
#############################################
##### compare between gender without locations as covariate #####
### number locations ####
loc_counts_female <- c(favstats(Choma5_clusts_shorter2_no_home_female$clust~Choma5_clusts_shorter2_no_home_female$partid)$n,
                       favstats(Mutasa5_clusts_shorter2_no_home_female$clust~Mutasa5_clusts_shorter2_no_home_female$partid)$n,
                       favstats(Nchelenge5_clusts_shorter2_no_home_female$clust~Nchelenge5_clusts_shorter2_no_home_female$partid)$n)
loc_counts_male <- c(favstats(Choma5_clusts_shorter2_no_home_male$clust~Choma5_clusts_shorter2_no_home_male$partid)$n,
                     favstats(Mutasa5_clusts_shorter2_no_home_male$clust~Mutasa5_clusts_shorter2_no_home_male$partid)$n,
                     favstats(Nchelenge5_clusts_shorter2_no_home_male$clust~Nchelenge5_clusts_shorter2_no_home_male$partid)$n)
only_home_female1 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_female$partid))
only_home_male1 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_male$partid))
only_home_female2 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_female$partid))
only_home_male2 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_male$partid))
only_home_female3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_female$partid))
only_home_male3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_male$partid))

loc_counts_female2 <- c(loc_counts_female,rep(0, only_home_female1+only_home_female2+only_home_female3))
loc_counts_male2 <- c(loc_counts_male,rep(0, only_home_male1 + only_home_male2 + only_home_male3))


locs_count2 <- as.data.frame(cbind(c(loc_counts_female2, loc_counts_male2),
                                   c(rep("female",length(loc_counts_female2)),
                                     rep("male",length(loc_counts_male2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-2,27)) +
  scale_x_continuous(breaks=c(seq(0,25,5)))

shapiro.test(locs_count2$values) ## not normal
kruskal.test(locs_count2$values, locs_count2$type) # p = 0.004

my_comp <- list(c("female","male"))
ggplot(data=locs_count2, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("red","blue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 30) 

##
### distance ####
medians_km_female <- c((favstats(Choma5_clusts_shorter2_no_home_female$hhdist_m_haversine~Choma5_clusts_shorter2_no_home_female$partid)$median)/1000,
                       (favstats(Mutasa5_clusts_shorter2_no_home_female$hhdist_m_haversine~Mutasa5_clusts_shorter2_no_home_female$partid)$median)/1000, 
                       (favstats(Nchelenge5_clusts_shorter2_no_home_female$hhdist_m_haversine2~Nchelenge5_clusts_shorter2_no_home_female$partid)$median)/1000)
medians_km_male <- c((favstats(Choma5_clusts_shorter2_no_home_male$hhdist_m_haversine~Choma5_clusts_shorter2_no_home_male$partid)$median)/1000,
                     (favstats(Mutasa5_clusts_shorter2_no_home_male$hhdist_m_haversine~Mutasa5_clusts_shorter2_no_home_male$partid)$median)/1000, 
                     (favstats(Nchelenge5_clusts_shorter2_no_home_male$hhdist_m_haversine2~Nchelenge5_clusts_shorter2_no_home_male$partid)$median)/1000)
favstats(medians_km_female) # km
favstats(medians_km_male) # km
medians_km_female_short <- medians_km_female[-which(medians_km_female > 50)]
favstats(medians_km_female_short)
medians_km_male_short <- medians_km_male[-which(medians_km_male > 50)]
favstats(medians_km_male_short)

dists_km2_gender <- as.data.frame(cbind(c(medians_km_female_short, medians_km_male_short),
                                        c(rep("female",length(medians_km_female_short)),
                                          rep("male",length(medians_km_male_short)))))
colnames(dists_km2_gender) <- c("values","type")
dists_km2_gender$values <- as.numeric(as.character(dists_km2_gender$values))
dists_km2_gender$type <- factor(dists_km2_gender$type)

ggplot(data=dists_km2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,43))


shapiro.test(dists_km2_gender$values) ## not normal
kruskal.test(dists_km2_gender$values, dists_km2_gender$type) # p = 0.6

my_comp <- list(c("female","male"))
ggplot(data=dists_km2_gender, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("red","blue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 48) 
#
### number of trips (without home) ####
trips_part_median_female <- c(favstats(Choma5_clust_trips_shorter2_no_home3_female$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_female$partid)$median,
                              favstats(Mutasa5_clust_trips_shorter2_no_home3_female$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_female$partid)$median,
                              favstats(Nchelenge5_clust_trips_shorter2_no_home3_female$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_female$partid)$median)
trips_part_median_male <- c(favstats(Choma5_clust_trips_shorter2_no_home3_male$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_male$partid)$median,
                            favstats(Mutasa5_clust_trips_shorter2_no_home3_male$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_male$partid)$median,
                            favstats(Nchelenge5_clust_trips_shorter2_no_home3_male$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_male$partid)$median)
favstats(trips_part_median_female)
#   min  Q1   median   Q3      max       mean       sd  n missing
#    1  1    1.5  2  13 1.654321 1.200781 162       0
favstats(trips_part_median_male)
#   min Q1 median  Q3   max      mean       sd  n missing
#    1  1      2  2  13 1.895105 1.470189 143       0

trips_part2_gender <- as.data.frame(cbind(c(trips_part_median_female, trips_part_median_male),
                                          c(rep("female",length(trips_part_median_female)),
                                            rep("male",length(trips_part_median_male)))))
colnames(trips_part2_gender) <- c("values","type")
trips_part2_gender$values <- as.numeric(as.character(trips_part2_gender$values))
trips_part2_gender$type <- factor(trips_part2_gender$type)

ggplot(data=trips_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,14)) + 
  scale_x_continuous(breaks=c(seq(0,14,2)))


shapiro.test(trips_part2_gender$values) ## not normal
kruskal.test(trips_part2_gender$values, trips_part2_gender$type) # p = 0.03

my_comp <- list(c("female","male"))
ggplot(data=trips_part2_gender, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("red","blue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 15) 
#
### time per trip ####
trip_time_part_median_female <- c(favstats(Choma5_clust_trips_shorter2_no_home_female$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_female$partid)$median,
                                  favstats(Mutasa5_clust_trips_shorter2_no_home_female$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_female$partid)$median,
                                  favstats(Nchelenge5_clust_trips_shorter2_no_home_female$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_female$partid)$median)
trip_time_part_median_male <- c(favstats(Choma5_clust_trips_shorter2_no_home_male$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_male$partid)$median,
                                favstats(Mutasa5_clust_trips_shorter2_no_home_male$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_male$partid)$median,
                                favstats(Nchelenge5_clust_trips_shorter2_no_home_male$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_male$partid)$median)
favstats(trip_time_part_median_female)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4844444 1.482014 2.30875 3.359931 33.96292 2.926667 3.12571 162       0
favstats(trip_time_part_median_male)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.4955556 1.319444 1.918194 2.885486 1218.114 10.96917 101.6743 143       0
trip_time_part_median_female2 <- trip_time_part_median_female[which(trip_time_part_median_female < 30)]
favstats(trip_time_part_median_female2)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.4844444 1.477083 2.303333 3.359722 12.22667 2.733895 1.942531 161       0
trip_time_part_median_male2 <- trip_time_part_median_male[which(trip_time_part_median_male < 30)]
favstats(trip_time_part_median_male2)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.4955556 1.315833 1.916181 2.876146 9.897222 2.468155 1.859071 142       0


trip_times_part2_gender <- as.data.frame(cbind(c(trip_time_part_median_female2, trip_time_part_median_male2),
                                               c(rep("female",length(trip_time_part_median_female2)),
                                                 rep("male",length(trip_time_part_median_male2)))))
colnames(trip_times_part2_gender) <- c("values","type")
trip_times_part2_gender$values <- as.numeric(as.character(trip_times_part2_gender$values))
trip_times_part2_gender$type <- factor(trip_times_part2_gender$type)

ggplot(data=trip_times_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,12)) +
  scale_x_continuous(breaks=c(seq(0,14,2)))
#

shapiro.test(trip_times_part2_gender$values) ## not normal
kruskal.test(trip_times_part2_gender$values, trip_times_part2_gender$type) # p = 0.05

my_comp <- list(c("female","male"))
ggplot(data=trip_times_part2_gender, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("red","blue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 15) 
#
#
### percent time at home ####
perc_female <- c(Choma5_clusts_shorter2_home_female$perc_time_home,
                 Mutasa5_clusts_shorter2_home_female$perc_time_home,
                 Nchelenge5_clusts_shorter2_home_female$perc_time_home)
perc_male <- c(Choma5_clusts_shorter2_home_male$perc_time_home,
                 Mutasa5_clusts_shorter2_home_male$perc_time_home,
                 Nchelenge5_clusts_shorter2_home_male$perc_time_home)

perc_time_home2_gender <- as.data.frame(cbind(c(perc_female, perc_male),
                                              c(rep("female",length(perc_female)),
                                                rep("male",length(perc_male)))))
colnames(perc_time_home2_gender) <- c("values","type")
perc_time_home2_gender$values <- as.numeric(as.character(perc_time_home2_gender$values))
perc_time_home2_gender$type <- factor(perc_time_home2_gender$type)

ggplot(data=perc_time_home2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,100))


shapiro.test(perc_time_home2_gender$values) ## not normal
kruskal.test(perc_time_home2_gender$values, perc_time_home2_gender$type) # p = 0.02

my_comp <- list(c("female","male"))
ggplot(data=perc_time_home2_gender, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("red","blue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 115) 
#
#
### Average Direction of Travel #####
Cho_ind <- which(colnames(Choma5_clusts_shorter2_no_home) %in% c("partid", "male", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))
Mut_ind <- which(colnames(Mutasa5_clusts_shorter2_no_home) %in% c("partid", "male", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))
Nch_ind <- which(colnames(Nchelenge5_clusts_shorter2_no_home) %in% c("partid", "male", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))

temp <- rbind(Choma5_clusts_shorter2_no_home[!duplicated(Choma5_clusts_shorter2_no_home$partid),c(Cho_ind)],
          Mutasa5_clusts_shorter2_no_home[!duplicated(Mutasa5_clusts_shorter2_no_home$partid),c(Mut_ind)],
          Nchelenge5_clusts_shorter2_no_home[!duplicated(Nchelenge5_clusts_shorter2_no_home$partid),c(Nch_ind)])
temp2 <- temp[which(!(is.na(temp$male))),]

favstats(temp2$avg_dir_from_home_bearing~temp2$male)
# temp2$male      min       Q1    median       Q3      max     mean        sd  n missing
# 1          0 2.0642705 96.63572 195.5242 287.6080 358.9609 190.5493 101.8336 162       0
# 2          1 0.3598077 93.74533 205.8879 275.4374 359.7089 188.0893 101.6469 143       0
favstats(temp2$avg_dir_from_home_mag~temp2$male)
# temp2$male       min        Q1    median       Q3      max     mean       sd  n missing
# 1          0 0.03380367 0.6581526 1.699842 9.977547 210.8254 10.41685 27.52569 162       0
# 2          1 0.09139251 0.6843518 1.556306 8.305963 205.0693 11.71859 29.46810 143       0

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=male), position=position_dodge2(preserve = "single"), width=2) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Male", values=c("red","blue"))


temp3 <- temp2[which(temp2$avg_dir_from_home_mag < 100),] ## 164 of 176 parts
favstats(temp3$avg_dir_from_home_bearing~temp3$male)
# temp3$male      min       Q1    median       Q3      max     mean        sd  n missing
# 1          0 2.0642705 92.43699 193.1252 285.8122 358.9609 188.1360 102.2103 157       0
# 2          1 0.3598077 94.71717 204.8540 276.3744 359.7089 188.7877 101.6721 140       0
favstats(temp3$avg_dir_from_home_mag~temp3$male)
# temp3$male       min        Q1    median       Q3      max     mean       sd  n missing
# 1          0 0.03380367 0.6530909 1.663715 7.658468 53.84515 6.034424  9.742668 157       0
# 2          1 0.09139251 0.6515056 1.495520 6.315138 87.08940 8.216004 16.832956 140       0

ggplot(data=temp3, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=male), position=position_dodge2(preserve = "single"), width=1) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Male", values=c("red","blue"))


ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=..count.., fill=male)) +
  geom_histogram(position = "dodge", binwidth = 22.5, alpha=0.8) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~male, labeller = labeller(male = c("0"="Female", "1"="Male")))+
  scale_fill_manual(values=c("red","blue"))


temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)] <- -(360-temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)])
ggplot(data=temp2, aes(x=avg_dir_from_home_values, y=..count.., fill=male)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~male, labeller = labeller(male = c("0"="Female", "1"="Male")))+
  scale_fill_manual(values=c("red","blue"))


temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values <0)] <- (360+temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values < 0)])

shapiro.test(temp2$avg_dir_from_home_values) ## not normal
kruskal.test(temp2$avg_dir_from_home_values, temp2$male) # p = 0.08

my_comp <- list(c("0","1"))
ggplot(data=temp2, aes(x=male, y=avg_dir_from_home_values)) +
  geom_boxplot(aes(fill=male)) +
  theme_classic() +
  scale_fill_manual(values = c("red","blue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 400) 






## average direction calc
temp2$values2 <- 90 - temp2$avg_dir_from_home_values
temp2$values2[which(temp2$values2 < 0)] <- temp2$values2[which(temp2$values2 < 0)] + 360
temp2$values3 <- deg2rad(temp2$values2)

temp2$avg_value <- as.numeric(NA)
for(ii in 1:nlevels(temp2$male)){
  level <- levels(temp2$male)[ii]
  temp <- temp2[which(temp2$male == level),]
  sum_cos_dist <- sum(cos(temp$values3))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$values3))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  temp2$avg_value[which(temp2$male == level)] <- avg_deg
}

temp2$avg_value2 <- 90 - temp2$avg_value
temp2$avg_value2[which(temp2$avg_value2 < 0)] <- temp2$avg_value2[which(temp2$avg_value2 < 0)] + 360
favstats(temp2$avg_value2~temp2$male)
## female: 338.52 --> NNW
## male: 290.95 --> WNW








#
#
### Time at locations ####
clust_time_part_median_female <- c((favstats(Choma5_clusts_shorter2_no_home_female$clust_trips_time_new~Choma5_clusts_shorter2_no_home_female$partid)$median)*24,
                                   (favstats(Mutasa5_clusts_shorter2_no_home_female$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_female$partid)$median)*24,
                                   (favstats(Nchelenge5_clusts_shorter2_no_home_female$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_female$partid)$median)*24)
clust_time_part_median_male <- c((favstats(Choma5_clusts_shorter2_no_home_male$clust_trips_time_new~Choma5_clusts_shorter2_no_home_male$partid)$median)*24,
                                 (favstats(Mutasa5_clusts_shorter2_no_home_male$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_male$partid)$median)*24,
                                 (favstats(Nchelenge5_clusts_shorter2_no_home_male$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_male$partid)$median)*24)
favstats(clust_time_part_median_female)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4605556 2.179931 3.267153 4.693403 98.0725 4.987047 9.696488 162       0
favstats(clust_time_part_median_male)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.7801389 2.054375   3.27 4.769375 1218.114 15.50819 103.2215 143       0
clust_time_part_median_female2 <- clust_time_part_median_female[which(clust_time_part_median_female < 35)]
clust_time_part_median_male2 <- clust_time_part_median_male[which(clust_time_part_median_male < 35)]
favstats(clust_time_part_median_female2)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4605556 2.111806 3.240139 4.596181 19.32806 3.761438 2.456596 159       0
favstats(clust_time_part_median_male2)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.7801389 2.048854 3.231944 4.549271 22.93264 4.109886 3.490424 138       0

clust_time_part_gender <- as.data.frame(cbind(c(clust_time_part_median_female2, clust_time_part_median_male2),
                                              c(rep("female",length(clust_time_part_median_female2)),
                                                rep("male",length(clust_time_part_median_male2)))))
colnames(clust_time_part_gender) <- c("values","type")
clust_time_part_gender$values <- as.numeric(as.character(clust_time_part_gender$values))
clust_time_part_gender$type <- factor(clust_time_part_gender$type)

ggplot(data=clust_time_part_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,24)) +
  scale_x_continuous(breaks=c(seq(0,30,2)))
#

shapiro.test(clust_time_part_gender$values) ## not normal
kruskal.test(clust_time_part_gender$values, clust_time_part_gender$type) # p = 1

my_comp <- list(c("female","male"))
ggplot(data=clust_time_part_gender, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("red","blue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 25) 
#
#
#### Biting time metrics #####
### number locations ####
loc_counts_female <- c(favstats(Choma5_clusts_shorter2_no_home_biting_places2_female$clust~Choma5_clusts_shorter2_no_home_biting_places2_female$partid)$n,
                       favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_female$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_female$partid)$n,
                       favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_female$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_female$partid)$n)
loc_counts_male <- c(favstats(Choma5_clusts_shorter2_no_home_biting_places2_male$clust~Choma5_clusts_shorter2_no_home_biting_places2_male$partid)$n,
                       favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_male$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_male$partid)$n,
                       favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_male$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_male$partid)$n)
only_home_female1 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_female$partid))
only_home_male1 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_male$partid))
only_home_female2 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_female$partid))
only_home_male2 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_male$partid))
only_home_female3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$male == 0)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_female$partid))
only_home_male3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$male == 1)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_male$partid))

loc_counts_female2 <- c(loc_counts_female,rep(0, only_home_female1+only_home_female2+only_home_female3))
loc_counts_male2 <- c(loc_counts_male,rep(0, only_home_male1 + only_home_male2 + only_home_male3))

locs_count2 <- as.data.frame(cbind(c(loc_counts_female2, loc_counts_male2),
                                   c(rep("female",length(loc_counts_female2)),
                                     rep("male",length(loc_counts_male2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(-0.5,13)) +
  scale_x_continuous(breaks=c(seq(0,25,2)))
##
### time per location ####
clust_time_part_median_female <- c((favstats(Choma5_clusts_shorter2_no_home_biting_places2_female$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_female$partid)$median)*24,
                                   (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_female$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_female$partid)$median)*24,
                                   (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_female$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_female$partid)$median)*24)
clust_time_part_median_male <- c((favstats(Choma5_clusts_shorter2_no_home_biting_places2_male$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_male$partid)$median)*24,
                                 (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_male$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_male$partid)$median)*24,
                                 (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_male$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_male$partid)$median)*24)
clust_time_part_median_male2 <- clust_time_part_median_male[which(clust_time_part_median_male < 500)]

clust_time_part_gender <- as.data.frame(cbind(c(clust_time_part_median_female, clust_time_part_median_male2),
                                              c(rep("female",length(clust_time_part_median_female)),
                                                rep("male",length(clust_time_part_median_male2)))))
colnames(clust_time_part_gender) <- c("values","type")
clust_time_part_gender$values <- as.numeric(as.character(clust_time_part_gender$values))
clust_time_part_gender$type <- factor(clust_time_part_gender$type)

ggplot(data=clust_time_part_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=4) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,300)) +
  scale_x_continuous(breaks=c(seq(0,300,50)))
##
### number trips ####
### PER PART ###
trips_part_median_female <- c(favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_female$partid)$median,
                              favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_female$partid)$median,
                              favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_female$partid)$median)
trips_part_median_male <- c(favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_male$partid)$median,
                              favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_male$partid)$median,
                              favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_male$partid)$median)
trips_part_median_male2 <- trips_part_median_male[which(trips_part_median_male < 20)]


trips_part2_gender <- as.data.frame(cbind(c(trips_part_median_female, trips_part_median_male2),
                                          c(rep("female",length(trips_part_median_female)),
                                            rep("male",length(trips_part_median_male2)))))
colnames(trips_part2_gender) <- c("values","type")
trips_part2_gender$values <- as.numeric(as.character(trips_part2_gender$values))
trips_part2_gender$type <- factor(trips_part2_gender$type)

ggplot(data=trips_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,14)) + 
  scale_x_continuous(breaks=c(seq(0,14,2)))



female_all <- c(Choma5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2, 
                Mutasa5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2,
                Nchelenge5_clust_trips_shorter2_no_home_biting_places3_female$trip_count_new2)
male_all <- c(Choma5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2, 
                Mutasa5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2,
                Nchelenge5_clust_trips_shorter2_no_home_biting_places3_male$trip_count_new2)

trips_all2 <- as.data.frame(cbind(c(female_all, male_all),
                                  c(rep("female",length(female_all)),
                                    rep("male",length(male_all)))))
colnames(trips_all2) <- c("values","type")
trips_all2$values <- as.numeric(as.character(trips_all2$values))
trips_all2$type <- factor(trips_all2$type)

table(trips_all2$type, as.factor(trips_all2$values))
#
##### percent biting time spent at home vs elsewhere #####
perc_female <- c(Choma5_clusts_shorter2_biting_places3b_female$percent_home_biting_time, 
                 Mutasa5_clusts_shorter2_biting_places3b_female$percent_home_biting_time, 
                 Nchelenge5_clusts_shorter2_biting_places3b_female$percent_home_biting_time)
perc_male <- c(Choma5_clusts_shorter2_biting_places3b_male$percent_home_biting_time, 
               Mutasa5_clusts_shorter2_biting_places3b_male$percent_home_biting_time, 
               Nchelenge5_clusts_shorter2_biting_places3b_male$percent_home_biting_time)

perc_time_home2_gender <- as.data.frame(cbind(c(perc_female, perc_male),
                                              c(rep("female",length(perc_female)),
                                                rep("male",length(perc_male)))))
colnames(perc_time_home2_gender) <- c("values","type")
perc_time_home2_gender$values <- as.numeric(as.character(perc_time_home2_gender$values))
perc_time_home2_gender$type <- factor(perc_time_home2_gender$type)

ggplot(data=perc_time_home2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(xlim=c(0,100))
##
#############################################
#############################################
##### compare between age without locations as covariate #####
### number of locations (without home) ####
loc_counts_1_16 <- c(favstats(Choma5_clusts_shorter2_no_home_1_16$clust~Choma5_clusts_shorter2_no_home_1_16$partid)$n,
                     favstats(Mutasa5_clusts_shorter2_no_home_1_16$clust~Mutasa5_clusts_shorter2_no_home_1_16$partid)$n,
                     favstats(Nchelenge5_clusts_shorter2_no_home_1_16$clust~Nchelenge5_clusts_shorter2_no_home_1_16$partid)$n)
loc_counts_17_34 <- c(favstats(Choma5_clusts_shorter2_no_home_17_34$clust~Choma5_clusts_shorter2_no_home_17_34$partid)$n,
                      favstats(Mutasa5_clusts_shorter2_no_home_17_34$clust~Mutasa5_clusts_shorter2_no_home_17_34$partid)$n,
                      favstats(Nchelenge5_clusts_shorter2_no_home_17_34$clust~Nchelenge5_clusts_shorter2_no_home_17_34$partid)$n)
loc_counts_35_54 <- c(favstats(Choma5_clusts_shorter2_no_home_35_54$clust~Choma5_clusts_shorter2_no_home_35_54$partid)$n,
                      favstats(Mutasa5_clusts_shorter2_no_home_35_54$clust~Mutasa5_clusts_shorter2_no_home_35_54$partid)$n,
                      favstats(Nchelenge5_clusts_shorter2_no_home_35_54$clust~Nchelenge5_clusts_shorter2_no_home_35_54$partid)$n)
loc_counts_55_up <- c(favstats(Choma5_clusts_shorter2_no_home_55_up$clust~Choma5_clusts_shorter2_no_home_55_up$partid)$n,
                      favstats(Mutasa5_clusts_shorter2_no_home_55_up$clust~Mutasa5_clusts_shorter2_no_home_55_up$partid)$n,
                      favstats(Nchelenge5_clusts_shorter2_no_home_55_up$clust~Nchelenge5_clusts_shorter2_no_home_55_up$partid)$n)

only_home_1_161 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_1_16)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_1_16$partid))
only_home_17_341 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_17_34)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_17_34$partid))
only_home_35_541 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_35_54)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_35_54$partid))
only_home_55_up1 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_55_up)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_55_up$partid))
only_home_1_162 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_1_16)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_1_16$partid))
only_home_17_342 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_17_34)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_17_34$partid))
only_home_35_542 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_35_54)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_35_54$partid))
only_home_55_up2 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_55_up)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_55_up$partid))
only_home_1_163 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_1_16)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_1_16$partid))
only_home_17_343 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_17_34)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_17_34$partid))
only_home_35_543 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_35_54)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_35_54$partid))
only_home_55_up3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_55_up)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_55_up$partid))

loc_counts_1_162 <- c(loc_counts_1_16,rep(0, only_home_1_161 + only_home_1_162 + only_home_1_163))
loc_counts_17_342 <- c(loc_counts_17_34,rep(0, only_home_17_341 + only_home_17_342 + only_home_17_343))
loc_counts_35_542 <- c(loc_counts_35_54,rep(0, only_home_35_541 + only_home_35_542 + only_home_35_543))
loc_counts_55_up2 <- c(loc_counts_55_up,rep(0, only_home_55_up1 + only_home_55_up2 + only_home_55_up3))

favstats(loc_counts_1_162)
# min   Q1 median     Q3  max     mean        sd  n missing
#   0  3      8 9.5  18 7.511628 4.973128 43       0
favstats(loc_counts_17_342)
# min   Q1 median     Q3  max     mean        sd  n missing
#    0  4      7 12.25  25 8.625 6.155822 112       0
favstats(loc_counts_35_542)
# min   Q1 median     Q3  max     mean        sd  n missing
#  0  5      9 14  22 9.439252 5.823144 107       0
favstats(loc_counts_55_up2)
# min   Q1 median     Q3  max     mean        sd  n missing
#   0  5      9 13  24 9.22807 6.032701 57       0

locs_count2 <- as.data.frame(cbind(c(loc_counts_1_162, loc_counts_17_342, loc_counts_35_542, loc_counts_55_up2),
                                   c(rep("1-16",length(loc_counts_1_162)),
                                     rep("17-34",length(loc_counts_17_342)),
                                     rep("35-54",length(loc_counts_35_542)),
                                     rep("55+",length(loc_counts_55_up2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-2,27)) +
  scale_x_continuous(breaks=c(seq(0,25,5)))


shapiro.test(locs_count2$values) ## not normal
kruskal.test(locs_count2$values, locs_count2$type) # p = 0.3

my_comp <- list(c("1-16", "17-34"), c("1-16","35-54"), c("1-16","55+"), c("17-34","35-54"), c("17-34","55+"), c("35-54","55+"))
ggplot(data=locs_count2, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 42) 

#
### distance of locations (without home) ####
medians_km_1_16 <- c((favstats(Choma5_clusts_shorter2_no_home_1_16$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_1_16$partid)$median),
                     (favstats(Mutasa5_clusts_shorter2_no_home_1_16$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_1_16$partid)$median),
                     (favstats(Nchelenge5_clusts_shorter2_no_home_1_16$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_1_16$partid)$median))
medians_km_17_34 <- c((favstats(Choma5_clusts_shorter2_no_home_17_34$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_17_34$partid)$median),
                      (favstats(Mutasa5_clusts_shorter2_no_home_17_34$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_17_34$partid)$median),
                      (favstats(Nchelenge5_clusts_shorter2_no_home_17_34$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_17_34$partid)$median))
medians_km_35_54 <- c((favstats(Choma5_clusts_shorter2_no_home_35_54$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_35_54$partid)$median),
                      (favstats(Mutasa5_clusts_shorter2_no_home_35_54$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_35_54$partid)$median),
                      (favstats(Nchelenge5_clusts_shorter2_no_home_35_54$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_35_54$partid)$median))
medians_km_55_up <- c((favstats(Choma5_clusts_shorter2_no_home_55_up$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_55_up$partid)$median),
                      (favstats(Mutasa5_clusts_shorter2_no_home_55_up$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_55_up$partid)$median),
                      (favstats(Nchelenge5_clusts_shorter2_no_home_55_up$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_55_up$partid)$median))

favstats(medians_km_1_16) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.08078032 0.7485426 1.153074 2.196252 25.71711 2.864173 5.232606 42       0
favstats(medians_km_17_34) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#    0.1012507 0.6350982 1.253317 3.255688 217.5585 8.60211 30.086 108       0
favstats(medians_km_35_54) # km
#          min       Q1   median       Q3      max    mean       sd  n missing
#  0.1602008 0.8927782 1.693878 11.38005 299.3132 15.34231 42.68301 102       0
favstats(medians_km_55_up) # km
#          min         Q1   median       Q3      max     mean       sd  n missing
#  0.1013454 0.557168 1.051781 2.235168 22.49014 3.491519 5.628475 56       0

medians_km_17_342 <- medians_km_17_34[which(medians_km_17_34 < 100)]
medians_km_35_542 <- medians_km_35_54[which(medians_km_35_54 < 100)]
favstats(medians_km_17_342) # km
# min        Q1   median       Q3      max     mean       sd   n missing
# 0.1012507 0.6307958 1.207019 3.140324 85.42321 4.758477 10.79431 106       0
# medians_km_35_542 <- medians_km_35_54[which(medians_km_35_54 < 100)]
favstats(medians_km_35_542) # km
# min        Q1   median       Q3      max     mean       sd  n missing
# 0.1602008 0.8765423 1.574526 7.741897 57.42812 6.709606 10.86486 97       0

dists_km2_age <- as.data.frame(cbind(c(medians_km_1_16, medians_km_17_342, medians_km_35_542, medians_km_55_up),
                                     c(rep("1-16",length(medians_km_1_16)),
                                       rep("17-34",length(medians_km_17_342)),
                                       rep("35-54",length(medians_km_35_542)),
                                       rep("55+",length(medians_km_55_up)))))
colnames(dists_km2_age) <- c("values","type")
dists_km2_age$values <- as.numeric(as.character(dists_km2_age$values))
dists_km2_age$type <- factor(dists_km2_age$type)

ggplot(data=dists_km2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=10) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-2,87))

shapiro.test(dists_km2_age$values) ## not normal
kruskal.test(dists_km2_age$values, dists_km2_age$type) # p = 0.1
pairwise.wilcox.test(dists_km2_age$values, dists_km2_age$type, p.adjust.method = "holm") ## 
#         1-16 17-34 35-54
#   17-34 1.0  -     -    
#   35-54 0.4  0.4   -    
#   55+   1.0  1.0   0.2   

my_comp <- list(c("1-16", "17-34"), c("1-16","35-54"), c("1-16","55+"), c("17-34","35-54"), c("17-34","55+"), c("35-54","55+"))
ggplot(data=dists_km2_age, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 150) 

###
### number of trips (without home) ####
trips_part_median_1_16 <- c(favstats(Choma5_clust_trips_shorter2_no_home3_1_16$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_1_16$partid)$median,
                            favstats(Mutasa5_clust_trips_shorter2_no_home3_1_16$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_1_16$partid)$median,
                            favstats(Nchelenge5_clust_trips_shorter2_no_home3_1_16$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_1_16$partid)$median)
trips_part_median_17_34 <- c(favstats(Choma5_clust_trips_shorter2_no_home3_17_34$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_17_34$partid)$median,
                             favstats(Mutasa5_clust_trips_shorter2_no_home3_17_34$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_17_34$partid)$median,
                             favstats(Nchelenge5_clust_trips_shorter2_no_home3_17_34$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_17_34$partid)$median)
trips_part_median_35_54 <- c(favstats(Choma5_clust_trips_shorter2_no_home3_35_54$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_35_54$partid)$median,
                             favstats(Mutasa5_clust_trips_shorter2_no_home3_35_54$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_35_54$partid)$median,
                             favstats(Nchelenge5_clust_trips_shorter2_no_home3_35_54$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_35_54$partid)$median)
trips_part_median_55_up <- c(favstats(Choma5_clust_trips_shorter2_no_home3_55_up$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3_55_up$partid)$median,
                             favstats(Mutasa5_clust_trips_shorter2_no_home3_55_up$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3_55_up$partid)$median,
                             favstats(Nchelenge5_clust_trips_shorter2_no_home3_55_up$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_55_up$partid)$median)

favstats(trips_part_median_1_16)
#   min  Q1   median   Q3      max       mean       sd  n missing
#    1  1      2  2  13 2.071429 1.901769 42       0
favstats(trips_part_median_17_34)
#   min Q1 median  Q3   max      mean       sd  n missing
#  1  1    1.5  2  13 1.675926 1.284914 108       0
favstats(trips_part_median_35_54)
#   min  Q1   median   Q3      max       mean       sd  n missing
#   1  1    1.5  2   7 1.593137 0.8044751 102       0
favstats(trips_part_median_55_up)
#   min Q1 median  Q3   max      mean       sd  n missing
#     1  1      2  2 11.5 2.017857 1.609529 56       0

trips_part2_age <- as.data.frame(cbind(c(trips_part_median_1_16, trips_part_median_17_34, trips_part_median_35_54, trips_part_median_55_up),
                                       c(rep("1-16",length(trips_part_median_1_16)),
                                         rep("17-34",length(trips_part_median_17_34)),
                                         rep("35-54",length(trips_part_median_35_54)),
                                         rep("55+",length(trips_part_median_55_up)))))
colnames(trips_part2_age) <- c("values","type")
trips_part2_age$values <- as.numeric(as.character(trips_part2_age$values))
trips_part2_age$type <- factor(trips_part2_age$type)

ggplot(data=trips_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,14)) + 
  scale_x_continuous(breaks=c(seq(0,14,2)))


shapiro.test(trips_part2_age$values) ## not normal
kruskal.test(trips_part2_age$values, trips_part2_age$type) # p = 0.08
pairwise.wilcox.test(trips_part2_age$values, trips_part2_age$type, p.adjust.method = "holm") ## 
#         1-16 17-34 35-54
#   17-34 0.3  -     -    
#   35-54 0.3  1.0   -    
#   55+   1.0  0.3   0.3 

my_comp <- list(c("1-16", "17-34"), c("1-16","35-54"), c("1-16","55+"), c("17-34","35-54"), c("17-34","55+"), c("35-54","55+"))
ggplot(data=trips_part2_age, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 22) 

#
### time per trip ####
trip_time_part_median_1_16 <- c(favstats(Choma5_clust_trips_shorter2_no_home_1_16$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_1_16$partid)$median,
                                favstats(Mutasa5_clust_trips_shorter2_no_home_1_16$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_1_16$partid)$median,
                                favstats(Nchelenge5_clust_trips_shorter2_no_home_1_16$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_1_16$partid)$median)
trip_time_part_median_17_34 <- c(favstats(Choma5_clust_trips_shorter2_no_home_17_34$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_17_34$partid)$median,
                                 favstats(Mutasa5_clust_trips_shorter2_no_home_17_34$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_17_34$partid)$median,
                                 favstats(Nchelenge5_clust_trips_shorter2_no_home_17_34$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_17_34$partid)$median)
trip_time_part_median_35_54 <- c(favstats(Choma5_clust_trips_shorter2_no_home_35_54$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_35_54$partid)$median,
                                 favstats(Mutasa5_clust_trips_shorter2_no_home_35_54$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_35_54$partid)$median,
                                 favstats(Nchelenge5_clust_trips_shorter2_no_home_35_54$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_35_54$partid)$median)
trip_time_part_median_55_up <- c(favstats(Choma5_clust_trips_shorter2_no_home_55_up$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home_55_up$partid)$median,
                                 favstats(Mutasa5_clust_trips_shorter2_no_home_55_up$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home_55_up$partid)$median,
                                 favstats(Nchelenge5_clust_trips_shorter2_no_home_55_up$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_55_up$partid)$median)

favstats(trip_time_part_median_1_16)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4844444 1.414722 2.036319 3.697083 9.897222 2.95995 2.339644 42       0
favstats(trip_time_part_median_17_34)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.4966667 1.378819 2.062292 2.951458 1218.114 13.71346 116.9894 108       0
favstats(trip_time_part_median_35_54)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4955556 1.440035 2.148681 3.190347 33.96292 2.958272 3.726842 102       0
favstats(trip_time_part_median_55_up)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.6075 1.525486 2.146597 3.121944 7.917222 2.579169 1.485692 56       0

trip_time_part_median_17_34_short <- trip_time_part_median_17_34[which(trip_time_part_median_17_34 < 30)] ## no outliers to remove
favstats(trip_time_part_median_17_34_short)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.4966667 1.374722 2.039861 2.902986 12.22667 2.457386 1.738164 107       0
trip_time_part_median_35_54_short <- trip_time_part_median_35_54[which(trip_time_part_median_35_54 < 30)] ## no outliers to remove
favstats(trip_time_part_median_35_54_short)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0..4955556 1.44 2.137361 3.173889 12.13694 2.651295 2.078507 101       0


trip_times_part2_age <- as.data.frame(cbind(c(trip_time_part_median_1_16, trip_time_part_median_17_34_short, trip_time_part_median_35_54_short, trip_time_part_median_55_up),
                                            c(rep("1-16",length(trip_time_part_median_1_16)),
                                              rep("17-34",length(trip_time_part_median_17_34_short)),
                                              rep("35-54",length(trip_time_part_median_35_54_short)),
                                              rep("55+",length(trip_time_part_median_55_up)))))
colnames(trip_times_part2_age) <- c("values","type")
trip_times_part2_age$values <- as.numeric(as.character(trip_times_part2_age$values))
trip_times_part2_age$type <- factor(trip_times_part2_age$type)

ggplot(data=trip_times_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,12)) + 
  scale_x_continuous(breaks=c(seq(0,14,2)))


shapiro.test(trip_times_part2_age$values) ## not normal
kruskal.test(trip_times_part2_age$values, trip_times_part2_age$type) # p = 0.8
pairwise.wilcox.test(trip_times_part2_age$values, trip_times_part2_age$type, p.adjust.method = "holm") ## 
#         1-16 17-34 35-54
#   17-34 1  -     -    
#   35-54 1  1   -    
#   55+   1  1   1 

my_comp <- list(c("1-16", "17-34"), c("1-16","35-54"), c("1-16","55+"), c("17-34","35-54"), c("17-34","55+"), c("35-54","55+"))
ggplot(data=trip_times_part2_age, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 22) 

#
### percent time at home ####
perc_1_16 <- c(Choma5_clusts_shorter2_home_1_16$perc_time_home,
               Mutasa5_clusts_shorter2_home_1_16$perc_time_home,
               Nchelenge5_clusts_shorter2_home_1_16$perc_time_home)
perc_17_34 <- c(Choma5_clusts_shorter2_home_17_34$perc_time_home,
                Mutasa5_clusts_shorter2_home_17_34$perc_time_home,
                Nchelenge5_clusts_shorter2_home_17_34$perc_time_home)
perc_35_54 <- c(Choma5_clusts_shorter2_home_35_54$perc_time_home,
                Mutasa5_clusts_shorter2_home_35_54$perc_time_home,
                Nchelenge5_clusts_shorter2_home_35_54$perc_time_home)
perc_55_up <- c(Choma5_clusts_shorter2_home_55_up$perc_time_home,
                Mutasa5_clusts_shorter2_home_55_up$perc_time_home,
                Nchelenge5_clusts_shorter2_home_55_up$perc_time_home)

perc_home2_age <- as.data.frame(cbind(c(perc_1_16, perc_17_34, perc_35_54, perc_55_up),
                                      c(rep("1-16",length(perc_1_16)),
                                        rep("17-34",length(perc_17_34)),
                                        rep("35-54",length(perc_35_54)),
                                        rep("55+",length(perc_55_up)))))
colnames(perc_home2_age) <- c("values","type")
perc_home2_age$values <- as.numeric(as.character(perc_home2_age$values))
perc_home2_age$type <- factor(perc_home2_age$type)

ggplot(data=perc_home2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-1,102))

shapiro.test(perc_home2_age$values) ## not normal
kruskal.test(perc_home2_age$values, perc_home2_age$type) # p = 0.5
pairwise.wilcox.test(perc_home2_age$values, perc_home2_age$type, p.adjust.method = "holm") ## 
#         1-16 17-34 35-54
#   17-34 1  -     -    
#   35-54 1  1   -    
#   55+   0.7  1   1 

my_comp <- list(c("1-16", "17-34"), c("1-16","35-54"), c("1-16","55+"), c("17-34","35-54"), c("17-34","55+"), c("35-54","55+"))
ggplot(data=perc_home2_age, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 170) 

##
### Average Direction of Travel #####
Cho_ind <- which(colnames(Choma5_clusts_shorter2_no_home) %in% c("partid", "age_all", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))
Mut_ind <- which(colnames(Mutasa5_clusts_shorter2_no_home) %in% c("partid", "age_all", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))
Nch_ind <- which(colnames(Nchelenge5_clusts_shorter2_no_home) %in% c("partid", "age_all", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))

temp <- rbind(Choma5_clusts_shorter2_no_home[!duplicated(Choma5_clusts_shorter2_no_home$partid),c(Cho_ind)],
              Mutasa5_clusts_shorter2_no_home[!duplicated(Mutasa5_clusts_shorter2_no_home$partid),c(Mut_ind)],
              Nchelenge5_clusts_shorter2_no_home[!duplicated(Nchelenge5_clusts_shorter2_no_home$partid),c(Nch_ind)])
age_group <- temp$age_all[!duplicated(temp$partid)]

temp$age_group <- factor(0, levels=c(0:3)) # 0 = <16; 1 = 17-34; 2 = 35-54; 3= 55+
temp$age_group[which(age_group >= 17 & age_group < 35)] <- 1
temp$age_group[which(age_group >= 35 & age_group < 55)] <- 2
temp$age_group[which(age_group >= 55)] <- 3
summary(as.factor(temp$age_group))

favstats(temp$avg_dir_from_home_bearing~temp$age_group)
# temp$age_group       min        Q1   median       Q3      max     mean        sd   n missing
# 1              3 2.0642705  75.12332 149.5546 247.0076 341.6283 164.1259  97.80916  56       0
# 2              2 0.3598077 103.72631 211.3930 275.4477 358.9609 190.2461  98.50107 102       0
# 3              1 3.8889726 122.09999 202.8050 285.5527 359.7089 195.9071 100.27600 108       0
# 4              0 7.3662570  91.96050 226.2721 306.4442 357.0854 202.6975 114.34066  42       0
favstats(temp$avg_dir_from_home_mag~temp$age_group)
# temp$age_group        min        Q1   median        Q3       max      mean        sd   n missing
# 1              3 0.05659092 0.5553171 1.499280  4.746362  29.08639  4.440132  6.149066  56       0
# 2              2 0.07662000 0.7610804 2.587721 14.798775 205.06933 17.129233 36.500730 102       0
# 3              1 0.08965062 0.7719555 1.673646  8.250431 210.82539 11.462618 30.424610 108       0
# 4              0 0.03380367 0.4509354 1.178555  3.593856  22.15730  3.528334  5.508587  42       0

ggplot(data=temp, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=age_group), position=position_dodge2(preserve = "single"), width=4) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Age Group", values=agecols[c(5,4,2,1)], labels=c("55+", "35-54", "17-34", "0-16"))
#


temp2 <- temp[which(temp$avg_dir_from_home_mag < 100),]
favstats(temp2$avg_dir_from_home_bearing~temp2$age_group)
# temp2$age_group      min       Q1    median       Q3      max     mean        sd  n missing
# 1               3 2.0642705  75.12332 149.5546 247.0076 341.6283 164.1259  97.80916  56       0
# 2               2 0.3598077 101.89538 208.9850 273.8847 358.9609 189.4197  98.91335  96       0
# 3               1 3.8889726 121.58784 199.3766 282.0305 359.7089 194.1570 100.39675 106       0
# 4               0 7.3662570  91.96050 226.2721 306.4442 357.0854 202.6975 114.34066  42       0
favstats(temp2$avg_dir_from_home_mag~temp2$age_group)
# temp2$age_group       min        Q1    median        Q3      max      mean       sd  n
# 1               3 0.05659092 0.5553171 1.499280  4.746362 29.08639 4.440132  6.149066  56       0
# 2               2 0.07662000 0.7013618 1.840349 11.784475 87.08940 9.183348 15.729300  96       0
# 3               1 0.08965062 0.7545367 1.659634  6.487171 82.08658 7.904726 15.842823 106       0
# 4               0 0.03380367 0.4509354 1.178555  3.593856 22.15730 3.528334  5.508587  42       0

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=age_group), position=position_dodge2(preserve = "single"), width=2) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Age Group", values=agecols[c(5,4,2,1)], labels=c("55+", "35-54", "17-34", "0-16"))




temp$avg_dir_from_home_values[which(temp$avg_dir_from_home_values > 348.75)] <- -(360-temp$avg_dir_from_home_values[which(temp$avg_dir_from_home_values > 348.75)])
ggplot(data=temp, aes(x=avg_dir_from_home_values, y=..count.., fill=age_group)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~age_group, labeller = labeller(age_group = c("0"="0-16", "1"="17-34","2"="35-54", "3"="55+")))+
  scale_fill_manual(values=agecols[c(1,2,4,5)]) +
  scale_y_continuous(breaks=c(seq(0,12,2)))



temp$avg_dir_from_home_values[which(temp$avg_dir_from_home_values < 0)] <- (360+temp$avg_dir_from_home_values[which(temp$avg_dir_from_home_values < 0)])

shapiro.test(temp$avg_dir_from_home_values) ## not normal
kruskal.test(temp$avg_dir_from_home_values, temp$age_group) # p = 0.5
pairwise.wilcox.test(temp$avg_dir_from_home_values, temp$age_group, p.adjust.method = "holm") ## 
#         1-16 17-34 35-54
#   17-34 1  -     -    
#   35-54 1  1   -    
#   55+   1  0.8   1

my_comp <- list(c("0", "1"), c("0","2"), c("0","3"), c("1","2"), c("1","3"), c("2","3"))
ggplot(data=temp, aes(x=age_group, y=avg_dir_from_home_values)) +
  geom_boxplot(aes(fill=age_group)) +  
  theme_classic() +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 600) 












## average direction calc
temp2$values2 <- 90 - temp2$avg_dir_from_home_values
temp2$values2[which(temp2$values2 < 0)] <- temp2$values2[which(temp2$values2 < 0)] + 360
temp2$values3 <- deg2rad(temp2$values2)

temp2$avg_value <- as.numeric(NA)
for(ii in 1:nlevels(temp2$age_group)){
  level <- levels(temp2$age_group)[ii]
  temp <- temp2[which(temp2$age_group == level),]
  sum_cos_dist <- sum(cos(temp$values3))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$values3))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  temp2$avg_value[which(temp2$age_group == level)] <- avg_deg
}

temp2$avg_value2 <- 90 - temp2$avg_value
temp2$avg_value2[which(temp2$avg_value2 < 0)] <- temp2$avg_value2[which(temp2$avg_value2 < 0)] + 360
favstats(temp2$avg_value2~temp2$age_group)
## 1: 320.36 --> NW
## 2: 276.45 --> W
## 3: 283.49 --> WNW
## 4: 57.83 --> ENE

#
### Time at locations ####
clust_time_part_median_1_16 <- c((favstats(Choma5_clusts_shorter2_no_home_1_16$clust_trips_time_new~Choma5_clusts_shorter2_no_home_1_16$partid)$median)*24,
                                 (favstats(Mutasa5_clusts_shorter2_no_home_1_16$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_1_16$partid)$median)*24,
                                 (favstats(Nchelenge5_clusts_shorter2_no_home_1_16$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_1_16$partid)$median)*24)
clust_time_part_median_17_34 <- c((favstats(Choma5_clusts_shorter2_no_home_17_34$clust_trips_time_new~Choma5_clusts_shorter2_no_home_17_34$partid)$median)*24,
                                  (favstats(Mutasa5_clusts_shorter2_no_home_17_34$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_17_34$partid)$median)*24,
                                  (favstats(Nchelenge5_clusts_shorter2_no_home_17_34$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_17_34$partid)$median)*24)
clust_time_part_median_35_54 <- c((favstats(Choma5_clusts_shorter2_no_home_35_54$clust_trips_time_new~Choma5_clusts_shorter2_no_home_35_54$partid)$median)*24,
                                  (favstats(Mutasa5_clusts_shorter2_no_home_35_54$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_35_54$partid)$median)*24,
                                  (favstats(Nchelenge5_clusts_shorter2_no_home_35_54$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_35_54$partid)$median)*24)
clust_time_part_median_55_up <- c((favstats(Choma5_clusts_shorter2_no_home_55_up$clust_trips_time_new~Choma5_clusts_shorter2_no_home_55_up$partid)$median)*24,
                                  (favstats(Mutasa5_clusts_shorter2_no_home_55_up$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_55_up$partid)$median)*24,
                                  (favstats(Nchelenge5_clusts_shorter2_no_home_55_up$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_55_up$partid)$median)*24)

favstats(clust_time_part_median_1_16)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 1.309583 2.075938 2.981944 4.287743 110.4386 7.051819 17.47148 42       0
favstats(clust_time_part_median_17_34)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.7801389 2.28 3.237917 4.848924 1218.114 17.20903 118.1802 108       0
favstats(clust_time_part_median_35_54)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.8297222 1.836111 3.238333 4.549375 98.0725 5.324926 11.34782 102       0
favstats(clust_time_part_median_55_up)
#       min        Q1   median       Q3      max     mean       sd  n missing
#  0.4605556 2.301111 3.457917 5.392326 87.48472 6.002805 11.77366 56       0

clust_time_part_median_1_162 <- clust_time_part_median_1_16[which(clust_time_part_median_1_16 < 35)]
clust_time_part_median_17_342 <- clust_time_part_median_17_34[which(clust_time_part_median_17_34 < 35)]
clust_time_part_median_35_542 <- clust_time_part_median_35_54[which(clust_time_part_median_35_54 < 35)]
clust_time_part_median_55_up2 <- clust_time_part_median_55_up[which(clust_time_part_median_55_up < 35)]

favstats(clust_time_part_median_1_162)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 1.309583 2.028715 2.9425 4.152535 18.57042 3.666403 2.96205 40       0
favstats(clust_time_part_median_17_342)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.7801389 2.271667 3.232778 4.560139 12.90486 3.773959 2.300069 105       0
favstats(clust_time_part_median_35_542)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.8297222 1.818472 3.213472 4.469792 17.23194 3.821175 2.942026 100       0
favstats(clust_time_part_median_55_up2)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.4605556 2.298889 3.334444 5.109722 22.93264 4.521316 4.000074 55       0

clust_time_part_age <- as.data.frame(cbind(c(clust_time_part_median_1_162, clust_time_part_median_17_342,clust_time_part_median_35_542,clust_time_part_median_55_up2),
                                           c(rep("1-16",length(clust_time_part_median_1_162)),
                                             rep("17-34",length(clust_time_part_median_17_342)),
                                             rep("35-54",length(clust_time_part_median_35_542)),
                                             rep("55+",length(clust_time_part_median_55_up2)))))
colnames(clust_time_part_age) <- c("values","type")
clust_time_part_age$values <- as.numeric(as.character(clust_time_part_age$values))
clust_time_part_age$type <- factor(clust_time_part_age$type)

ggplot(data=clust_time_part_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,25)) +
  scale_x_continuous(breaks=c(seq(0,25,2)))
#


shapiro.test(clust_time_part_age$values) ## not normal
kruskal.test(clust_time_part_age$values, clust_time_part_age$type) # p = 0.6
pairwise.wilcox.test(clust_time_part_age$values, clust_time_part_age$type, p.adjust.method = "holm") ## 
#         1-16 17-34 35-54
#   17-34 1  -     -    
#   35-54 1  1   -    
#   55+   1  1   1 

my_comp <- list(c("1-16", "17-34"), c("1-16","35-54"), c("1-16","55+"), c("17-34","35-54"), c("17-34","55+"), c("35-54","55+"))
ggplot(data=clust_time_part_age, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 40) 

#
#
#### Biting time metrics #####
### number locations ####
loc_counts_1_16 <- c(favstats(Choma5_clusts_shorter2_no_home_biting_places2_1_16$clust~Choma5_clusts_shorter2_no_home_biting_places2_1_16$partid)$n,
                     favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_1_16$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_1_16$partid)$n,
                     favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16$partid)$n)
loc_counts_17_34 <- c(favstats(Choma5_clusts_shorter2_no_home_biting_places2_17_34$clust~Choma5_clusts_shorter2_no_home_biting_places2_17_34$partid)$n,
                      favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_17_34$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_17_34$partid)$n,
                      favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34$partid)$n)
loc_counts_35_54 <- c(favstats(Choma5_clusts_shorter2_no_home_biting_places2_35_54$clust~Choma5_clusts_shorter2_no_home_biting_places2_35_54$partid)$n,
                      favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_35_54$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_35_54$partid)$n,
                      favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54$partid)$n)
loc_counts_55_up <- c(favstats(Choma5_clusts_shorter2_no_home_biting_places2_55_up$clust~Choma5_clusts_shorter2_no_home_biting_places2_55_up$partid)$n,
                      favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_55_up$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_55_up$partid)$n,
                      favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up$partid)$n)

only_home_1_161 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_1_16)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_1_16$partid))
only_home_17_341 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_17_34)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_17_34$partid))
only_home_35_541 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_35_54)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_35_54$partid))
only_home_55_up1 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$partid %in% Choma_parts_55_up)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_55_up$partid))
only_home_1_162 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_1_16)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_1_16$partid))
only_home_17_342 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_17_34)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_17_34$partid))
only_home_35_542 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_35_54)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_35_54$partid))
only_home_55_up2 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$partid %in% Mutasa_parts_55_up)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_55_up$partid))
only_home_1_163 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_1_16)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16$partid))
only_home_17_343 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_17_34)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34$partid))
only_home_35_543 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_35_54)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54$partid))
only_home_55_up3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$partid %in% Nchelenge_parts_55_up)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up$partid))

loc_counts_1_162 <- c(loc_counts_1_16,rep(0, only_home_1_161 + only_home_1_162 + only_home_1_163))
loc_counts_17_342 <- c(loc_counts_17_34,rep(0, only_home_17_341 + only_home_17_342 + only_home_17_343))
loc_counts_35_542 <- c(loc_counts_35_54,rep(0, only_home_35_541 + only_home_35_542 + only_home_35_543))
loc_counts_55_up2 <- c(loc_counts_55_up,rep(0, only_home_55_up1 + only_home_55_up2 + only_home_55_up3))

locs_count2 <- as.data.frame(cbind(c(loc_counts_1_162, loc_counts_17_342, loc_counts_35_542, loc_counts_55_up2),
                                   c(rep("1-16",length(loc_counts_1_162)),
                                     rep("17-34",length(loc_counts_17_342)),
                                     rep("35-54",length(loc_counts_35_542)),
                                     rep("55+",length(loc_counts_55_up2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.4) + 
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-0.6,12)) +
  scale_x_continuous(breaks=c(seq(0,25,2)))
#
##
### time per location ####
clust_time_part_median_1_16 <- c((favstats(Choma5_clusts_shorter2_no_home_biting_places2_1_16$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_1_16$partid)$median)*24,
                                 (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_1_16$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_1_16$partid)$median)*24,
                                 (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_1_16$partid)$median)*24)
clust_time_part_median_17_34 <- c((favstats(Choma5_clusts_shorter2_no_home_biting_places2_17_34$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_17_34$partid)$median)*24,
                                  (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_17_34$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_17_34$partid)$median)*24,
                                  (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_17_34$partid)$median)*24)
clust_time_part_median_35_54 <- c((favstats(Choma5_clusts_shorter2_no_home_biting_places2_35_54$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_35_54$partid)$median)*24,
                                  (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_35_54$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_35_54$partid)$median)*24,
                                  (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_35_54$partid)$median)*24)
clust_time_part_median_55_up <- c((favstats(Choma5_clusts_shorter2_no_home_biting_places2_55_up$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_55_up$partid)$median)*24,
                                  (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_55_up$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_55_up$partid)$median)*24,
                                  (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_55_up$partid)$median)*24)

clust_time_part_median_17_342 <- clust_time_part_median_17_34[which(clust_time_part_median_17_34 < 500)]


clust_time_part_age <- as.data.frame(cbind(c(clust_time_part_median_1_16, clust_time_part_median_17_342,clust_time_part_median_35_54,clust_time_part_median_55_up),
                                           c(rep("1-16",length(clust_time_part_median_1_16)),
                                             rep("17-34",length(clust_time_part_median_17_342)),
                                             rep("35-54",length(clust_time_part_median_35_54)),
                                             rep("55+",length(clust_time_part_median_55_up)))))
colnames(clust_time_part_age) <- c("values","type")
clust_time_part_age$values <- as.numeric(as.character(clust_time_part_age$values))
clust_time_part_age$type <- factor(clust_time_part_age$type)

ggplot(data=clust_time_part_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 50, alpha=0.8) +
  geom_density(aes(color=type), adjust=7) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-14,300)) +
  scale_x_continuous(breaks=c(seq(0,300,50)))
#
##
### number trips ####
### PER PART ###
trips_part_median_1_16 <- c(favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2 ~ Choma5_clust_trips_shorter2_no_home_biting_places3_1_16$partid)$median,
                            favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2 ~ Mutasa5_clust_trips_shorter2_no_home_biting_places3_1_16$partid)$median,
                            favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_biting_places3_1_16$partid)$median)
trips_part_median_17_34 <- c(favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2 ~ Choma5_clust_trips_shorter2_no_home_biting_places3_17_34$partid)$median,
                             favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2 ~ Mutasa5_clust_trips_shorter2_no_home_biting_places3_17_34$partid)$median,
                             favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_biting_places3_17_34$partid)$median)
trips_part_median_35_54 <- c(favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2 ~ Choma5_clust_trips_shorter2_no_home_biting_places3_35_54$partid)$median,
                             favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2 ~ Mutasa5_clust_trips_shorter2_no_home_biting_places3_35_54$partid)$median,
                             favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_biting_places3_35_54$partid)$median)
trips_part_median_55_up <- c(favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2 ~ Choma5_clust_trips_shorter2_no_home_biting_places3_55_up$partid)$median,
                             favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2 ~ Mutasa5_clust_trips_shorter2_no_home_biting_places3_55_up$partid)$median,
                             favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_biting_places3_55_up$partid)$median)
trips_part_median_1_162 <- trips_part_median_1_16[which(trips_part_median_1_16 < 20)]
trips_part_median_35_542 <- trips_part_median_35_54[which(trips_part_median_35_54 < 20)]

trips_part2_age <- as.data.frame(cbind(c(trips_part_median_1_162, trips_part_median_17_34, trips_part_median_35_542, trips_part_median_55_up),
                                       c(rep("1-16",length(trips_part_median_1_162)),
                                         rep("17-34",length(trips_part_median_17_34)),
                                         rep("35-54",length(trips_part_median_35_542)),
                                         rep("55+",length(trips_part_median_55_up)))))
colnames(trips_part2_age) <- c("values","type")
trips_part2_age$values <- as.numeric(as.character(trips_part2_age$values))
trips_part2_age$type <- factor(trips_part2_age$type)

ggplot(data=trips_part2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(0,14)) + 
  scale_x_continuous(breaks=c(seq(0,14,2)))




all_1_16 <- c(Choma5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2, 
                Mutasa5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2,
                Nchelenge5_clust_trips_shorter2_no_home_biting_places3_1_16$trip_count_new2)
all_17_34 <- c(Choma5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2, 
              Mutasa5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2,
              Nchelenge5_clust_trips_shorter2_no_home_biting_places3_17_34$trip_count_new2)
all_35_54 <- c(Choma5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2, 
              Mutasa5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2,
              Nchelenge5_clust_trips_shorter2_no_home_biting_places3_35_54$trip_count_new2)
all_55_up <- c(Choma5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2, 
              Mutasa5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2,
              Nchelenge5_clust_trips_shorter2_no_home_biting_places3_55_up$trip_count_new2)


trips_all2 <- as.data.frame(cbind(c(all_1_16, all_17_34, all_35_54, all_55_up),
                                  c(rep("1-16",length(all_1_16)),
                                    rep("17-34",length(all_17_34)),
                                    rep("35-54",length(all_35_54)),
                                    rep("55+",length(all_55_up)))))
colnames(trips_all2) <- c("values","type")
trips_all2$values <- as.numeric(as.character(trips_all2$values))
trips_all2$type <- factor(trips_all2$type)

table(trips_all2$type, as.factor(trips_all2$values))
#
##### percent biting time spent at home vs elsewhere #####
perc_female <- c(Choma5_clusts_shorter2_biting_places3b_female$percent_home_biting_time, 
                 Mutasa5_clusts_shorter2_biting_places3b_female$percent_home_biting_time, 
                 Nchelenge5_clusts_shorter2_biting_places3b_female$percent_home_biting_time)
perc_male <- c(Choma5_clusts_shorter2_biting_places3b_male$percent_home_biting_time, 
               Mutasa5_clusts_shorter2_biting_places3b_male$percent_home_biting_time, 
               Nchelenge5_clusts_shorter2_biting_places3b_male$percent_home_biting_time)

perc_1_16 <- c(Choma5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time,
               Mutasa5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time,
               Nchelenge5_clusts_shorter2_biting_places3b_1_16$percent_home_biting_time)
perc_17_34 <- c(Choma5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time,
                Mutasa5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time,
                Nchelenge5_clusts_shorter2_biting_places3b_17_34$percent_home_biting_time)
perc_35_54 <- c(Choma5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time,
                Mutasa5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time,
                Nchelenge5_clusts_shorter2_biting_places3b_35_54$percent_home_biting_time)
perc_55_up <- c(Choma5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time,
                Mutasa5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time,
                Nchelenge5_clusts_shorter2_biting_places3b_55_up$percent_home_biting_time)

perc_home2_age <- as.data.frame(cbind(c(perc_1_16, perc_17_34, perc_35_54, perc_55_up),
                                      c(rep("1-16",length(perc_1_16)),
                                        rep("17-34",length(perc_17_34)),
                                        rep("35-54",length(perc_35_54)),
                                        rep("55+",length(perc_55_up)))))
colnames(perc_home2_age) <- c("values","type")
perc_home2_age$values <- as.numeric(as.character(perc_home2_age$values))
perc_home2_age$type <- factor(perc_home2_age$type)

ggplot(data=perc_home2_age, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) +
  scale_fill_manual(values = agecols[c(1,2,4,5)]) +
  scale_color_manual(values = agecols[c(1,2,4,5)]) +
  coord_cartesian(xlim=c(-1,102))
##
#############################################
#############################################
##### compare between rainy/dry without locations as covariate #####
### number locations ####
loc_counts_dry <- c(favstats(Choma5_rain_clusts_shorter2_no_home_dry$clust~Choma5_rain_clusts_shorter2_no_home_dry$partid)$n,
                       favstats(Mutasa5_rain_clusts_shorter2_no_home_dry$clust~Mutasa5_rain_clusts_shorter2_no_home_dry$partid)$n,
                       favstats(Nchelenge5_clusts_shorter2_no_home_dry$clust~Nchelenge5_clusts_shorter2_no_home_dry$partid)$n)
loc_counts_rainy <- c(favstats(Choma5_rain_clusts_shorter2_no_home_rainy$clust~Choma5_rain_clusts_shorter2_no_home_rainy$partid)$n,
                     favstats(Mutasa5_rain_clusts_shorter2_no_home_rainy$clust~Mutasa5_rain_clusts_shorter2_no_home_rainy$partid)$n,
                     favstats(Nchelenge5_clusts_shorter2_no_home_rainy$clust~Nchelenge5_clusts_shorter2_no_home_rainy$partid)$n)
only_home_dry1 <- nlevels(as.factor(Choma5_rain_clusts_shorter2$partid[which(Choma5_rain_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Choma5_rain_clusts_shorter2_no_home_dry$partid))
only_home_rainy1 <- nlevels(as.factor(Choma5_rain_clusts_shorter2$partid[which(Choma5_rain_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Choma5_rain_clusts_shorter2_no_home_rainy$partid))
only_home_dry2 <- nlevels(as.factor(Mutasa5_rain_clusts_shorter2$partid[which(Mutasa5_rain_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Mutasa5_rain_clusts_shorter2_no_home_dry$partid))
only_home_rainy2 <- nlevels(as.factor(Mutasa5_rain_clusts_shorter2$partid[which(Mutasa5_rain_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Mutasa5_rain_clusts_shorter2_no_home_rainy$partid))
only_home_dry3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_dry$partid))
only_home_rainy3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_rainy$partid))

loc_counts_dry2 <- c(loc_counts_dry,rep(0, only_home_dry1+only_home_dry2+only_home_dry3))
loc_counts_rainy2 <- c(loc_counts_rainy,rep(0, only_home_rainy1 + only_home_rainy2 + only_home_rainy3))


locs_count2 <- as.data.frame(cbind(c(loc_counts_dry2, loc_counts_rainy2),
                                   c(rep("dry",length(loc_counts_dry2)),
                                     rep("rainy",length(loc_counts_rainy2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 5, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-2,27)) +
  scale_x_continuous(breaks=c(seq(0,25,5)))

shapiro.test(locs_count2$values) ## not normal
kruskal.test(locs_count2$values, locs_count2$type) # p = 0.5

my_comp <- list(c("rainy","dry"))
ggplot(data=locs_count2, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 28) 

##
### distance ####
medians_km_dry <- c((favstats(Choma5_rain_clusts_shorter2_no_home_dry$hhdist_m_haversine~Choma5_rain_clusts_shorter2_no_home_dry$partid)$median)/1000,
                    (favstats(Mutasa5_rain_clusts_shorter2_no_home_dry$hhdist_m_haversine~Mutasa5_rain_clusts_shorter2_no_home_dry$partid)$median)/1000,
                    (favstats(Nchelenge5_clusts_shorter2_no_home_dry$hhdist_m_haversine~Nchelenge5_clusts_shorter2_no_home_dry$partid)$median)/1000)
medians_km_rainy <- c((favstats(Choma5_rain_clusts_shorter2_no_home_rainy$hhdist_m_haversine~Choma5_rain_clusts_shorter2_no_home_rainy$partid)$median)/1000,
                      (favstats(Mutasa5_rain_clusts_shorter2_no_home_rainy$hhdist_m_haversine~Mutasa5_rain_clusts_shorter2_no_home_rainy$partid)$median)/1000,
                      (favstats(Nchelenge5_clusts_shorter2_no_home_rainy$hhdist_m_haversine~Nchelenge5_clusts_shorter2_no_home_rainy$partid)$median)/1000)
favstats(medians_km_dry) # km
favstats(medians_km_rainy) # km
medians_km_dry_short <- medians_km_dry[-which(medians_km_dry > 50)]
favstats(medians_km_dry_short)
medians_km_rainy_short <- medians_km_rainy[-which(medians_km_rainy > 50)]
favstats(medians_km_rainy_short)

dists_km2_gender <- as.data.frame(cbind(c(medians_km_dry_short, medians_km_rainy_short),
                                        c(rep("dry",length(medians_km_dry_short)),
                                          rep("rainy",length(medians_km_rainy_short)))))
colnames(dists_km2_gender) <- c("values","type")
dists_km2_gender$values <- as.numeric(as.character(dists_km2_gender$values))
dists_km2_gender$type <- factor(dists_km2_gender$type)

ggplot(data=dists_km2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2.5) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,40))

shapiro.test(dists_km2_gender$values) ## not normal
kruskal.test(dists_km2_gender$values, dists_km2_gender$type) # p = 0.004

my_comp <- list(c("rainy","dry"))
ggplot(data=dists_km2_gender, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 42) 

##### Kruskall test for Nchelenge: p=0.12 ; Mutasa: p=0.38 ; Choma: p=0.002

#
### number of trips (without home) ####
trips_part_median_dry <- c(favstats(Choma5_rain_clust_trips_shorter2_no_home3_dry$trip_count_new ~ Choma5_rain_clust_trips_shorter2_no_home3_dry$partid)$median,
                              favstats(Mutasa5_rain_clust_trips_shorter2_no_home3_dry$trip_count_new ~ Mutasa5_rain_clust_trips_shorter2_no_home3_dry$partid)$median,
                              favstats(Nchelenge5_clust_trips_shorter2_no_home3_dry$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_dry$partid)$median)
trips_part_median_rainy <- c(favstats(Choma5_rain_clust_trips_shorter2_no_home3_rainy$trip_count_new ~ Choma5_rain_clust_trips_shorter2_no_home3_rainy$partid)$median,
                            favstats(Mutasa5_rain_clust_trips_shorter2_no_home3_rainy$trip_count_new ~ Mutasa5_rain_clust_trips_shorter2_no_home3_rainy$partid)$median,
                            favstats(Nchelenge5_clust_trips_shorter2_no_home3_rainy$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3_rainy$partid)$median)
favstats(trips_part_median_dry)
#   min  Q1   median   Q3      max       mean       sd  n missing
#    1  1    1.5  2  13 1.759398 1.265213 133       0
favstats(trips_part_median_rainy)
#   min Q1 median  Q3   max      mean       sd  n missing
#   1  1    1.5  2  13 1.809353 1.510106 139       0

trips_part2_gender <- as.data.frame(cbind(c(trips_part_median_dry, trips_part_median_rainy),
                                          c(rep("dry",length(trips_part_median_dry)),
                                            rep("rainy",length(trips_part_median_rainy)))))
colnames(trips_part2_gender) <- c("values","type")
trips_part2_gender$values <- as.numeric(as.character(trips_part2_gender$values))
trips_part2_gender$type <- factor(trips_part2_gender$type)

ggplot(data=trips_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,14)) + 
  scale_x_continuous(breaks=c(seq(0,14,2)))


shapiro.test(trips_part2_gender$values) ## not normal
kruskal.test(trips_part2_gender$values, trips_part2_gender$type) # p = 0.9

my_comp <- list(c("rainy","dry"))
ggplot(data=trips_part2_gender, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 14) 

#
### time per trip ####
trip_time_part_median_dry <- c(favstats(Choma5_rain_clust_trips_shorter2_no_home_dry$trip_time_new_hour ~ Choma5_rain_clust_trips_shorter2_no_home_dry$partid)$median,
                                  favstats(Mutasa5_rain_clust_trips_shorter2_no_home_dry$trip_time_new_hour ~ Mutasa5_rain_clust_trips_shorter2_no_home_dry$partid)$median,
                                  favstats(Nchelenge5_clust_trips_shorter2_no_home_dry$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_dry$partid)$median)
trip_time_part_median_rainy <- c(favstats(Choma5_rain_clust_trips_shorter2_no_home_rainy$trip_time_new_hour ~ Choma5_rain_clust_trips_shorter2_no_home_rainy$partid)$median,
                                favstats(Mutasa5_rain_clust_trips_shorter2_no_home_rainy$trip_time_new_hour ~ Mutasa5_rain_clust_trips_shorter2_no_home_rainy$partid)$median,
                                favstats(Nchelenge5_clust_trips_shorter2_no_home_rainy$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home_rainy$partid)$median)
favstats(trip_time_part_median_dry)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4955556 1.214722 2.063611 2.889167 33.96292 2.610676 3.165904 133       0
favstats(trip_time_part_median_rainy)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.4844444 1.387361 2.022778 3.416667 12.13694 2.732409 2.133464 139       0
trip_time_part_median_dry2 <- trip_time_part_median_dry[which(trip_time_part_median_dry < 30)]
favstats(trip_time_part_median_dry2)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4955556 1.214028 2.042083 2.882292 8.066389 2.373159 1.593449 132       0

trip_times_part2_gender <- as.data.frame(cbind(c(trip_time_part_median_dry2, trip_time_part_median_rainy),
                                               c(rep("dry",length(trip_time_part_median_dry2)),
                                                 rep("rainy",length(trip_time_part_median_rainy)))))
colnames(trip_times_part2_gender) <- c("values","type")
trip_times_part2_gender$values <- as.numeric(as.character(trip_times_part2_gender$values))
trip_times_part2_gender$type <- factor(trip_times_part2_gender$type)

ggplot(data=trip_times_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,12.25)) +
  scale_x_continuous(breaks=c(seq(0,14,2)))


shapiro.test(trip_times_part2_gender$values) ## not normal
kruskal.test(trip_times_part2_gender$values, trip_times_part2_gender$type) # p = 0.3

my_comp <- list(c("rainy","dry"))
ggplot(data=trip_times_part2_gender, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 14) 

#
#
### percent time at home ####
perc_dry <- c(Choma5_rain_clusts_shorter2_home_dry$perc_time_home,
                 Mutasa5_rain_clusts_shorter2_home_dry$perc_time_home,
                 Nchelenge5_clusts_shorter2_home_dry$perc_time_home)
perc_rainy <- c(Choma5_rain_clusts_shorter2_home_rainy$perc_time_home,
               Mutasa5_rain_clusts_shorter2_home_rainy$perc_time_home,
               Nchelenge5_clusts_shorter2_home_rainy$perc_time_home)

perc_time_home2_gender <- as.data.frame(cbind(c(perc_dry, perc_rainy),
                                              c(rep("dry",length(perc_dry)),
                                                rep("rainy",length(perc_rainy)))))
colnames(perc_time_home2_gender) <- c("values","type")
perc_time_home2_gender$values <- as.numeric(as.character(perc_time_home2_gender$values))
perc_time_home2_gender$type <- factor(perc_time_home2_gender$type)

ggplot(data=perc_time_home2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-2,102))


shapiro.test(perc_time_home2_gender$values) ## not normal
kruskal.test(perc_time_home2_gender$values, perc_time_home2_gender$type) # p = 0.3

my_comp <- list(c("rainy","dry"))
ggplot(data=perc_time_home2_gender, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 110) 

#
### Average Direction of Travel #####
Cho_ind <- which(colnames(Choma5_rain_clusts_shorter2_no_home) %in% c("partid", "rainy", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))
Mut_ind <- which(colnames(Mutasa5_rain_clusts_shorter2_no_home) %in% c("partid", "rainy", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))
Nch_ind <- which(colnames(Nchelenge5_clusts_shorter2_no_home) %in% c("partid", "rainy", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))

temp <- rbind(Choma5_rain_clusts_shorter2_no_home[!duplicated(Choma5_rain_clusts_shorter2_no_home$partid),c(Cho_ind)],
              Mutasa5_rain_clusts_shorter2_no_home[!duplicated(Mutasa5_rain_clusts_shorter2_no_home$partid),c(Mut_ind)],
              Nchelenge5_clusts_shorter2_no_home[!duplicated(Nchelenge5_clusts_shorter2_no_home$partid),c(Nch_ind)])
temp2 <- temp[which(!(is.na(temp$rainy))),]

favstats(temp2$avg_dir_from_home_bearing~temp2$rainy)
# temp2$rainy      min       Q1    median       Q3      max     mean        sd  n missing
# 1           0 0.3598077 122.35606 220.2642 288.8807 345.7862 202.1780 101.2205 133       0
# 2           1 3.8889726  85.98204 192.2631 262.5723 359.7089 179.9521 100.2945 139       0
favstats(temp2$avg_dir_from_home_mag~temp2$rainy)
# temp2$rainy       min        Q1    median       Q3      max     mean       sd  n missing
# 1           0 0.03380367 0.5699304 1.466226  4.602691 205.0693 11.746218 29.98367 133       0
# 2           1 0.08965062 0.7615483 2.003211 11.402396 189.2363  9.284823 23.49419 139       0 

ggplot(data=temp2, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=rainy), position=position_dodge2(preserve = "single"), width=2) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("rainy", values=c("sienna", "skyblue"))


temp3 <- temp2[which(temp2$avg_dir_from_home_mag < 100),] ## 164 of 176 parts
favstats(temp3$avg_dir_from_home_bearing~temp3$rainy)
# temp3$rainy      min       Q1    median       Q3      max     mean        sd  n missing
# 1           0 0.3598077 122.34616 221.0076 289.6122 345.7862 202.3208 101.3186 128       0
# 2           1 3.8889726  84.39213 189.7363 261.4512 359.7089 178.8367 100.5393 137       0
favstats(temp3$avg_dir_from_home_mag~temp3$rainy)
# temp3$rainy       min        Q1    median       Q3      max     mean       sd  n missing
# 1           0 0.03380367 0.5637756 1.439013  4.318618 87.08940 6.841985 15.42526 128       0
# 2           1 0.08965062 0.7606124 1.999662 11.039998 63.79277 6.731726 10.16681 137       0

ggplot(data=temp3, aes(x=avg_dir_from_home_bearing, y=avg_dir_from_home_mag)) +
  geom_col(aes(fill=rainy), position=position_dodge2(preserve = "single"), width=1) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("rainy", values=c("sienna", "skyblue"))



temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)] <- -(360-temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values > 348.75)])
ggplot(data=temp2, aes(x=avg_dir_from_home_values, y=..count.., fill=rainy)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~rainy, labeller = labeller(rainy = c("0"="Dry", "1"="Rainy")))+
  scale_fill_manual(values=c("sienna","skyblue")) 


temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values < 0)] <- (360+temp2$avg_dir_from_home_values[which(temp2$avg_dir_from_home_values < 0)])

shapiro.test(temp2$avg_dir_from_home_values) ## not normal
kruskal.test(temp2$avg_dir_from_home_values, temp2$rainy) # p = 0.03

my_comp <- list(c("0","1"))
ggplot(data=temp2, aes(x=rainy, y=avg_dir_from_home_values)) +
  geom_boxplot(aes(fill=rainy)) +  
  theme_classic() +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 405) 








## average direction calc
temp2$values2 <- 90 - temp2$avg_dir_from_home_values
temp2$values2[which(temp2$values2 < 0)] <- temp2$values2[which(temp2$values2 < 0)] + 360
temp2$values3 <- deg2rad(temp2$values2)

temp2$avg_value <- as.numeric(NA)
for(ii in 1:nlevels(temp2$rainy)){
  level <- levels(temp2$rainy)[ii]
  temp <- temp2[which(temp2$rainy == level),]
  sum_cos_dist <- sum(cos(temp$values3))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$values3))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  temp2$avg_value[which(temp2$rainy == level)] <- avg_deg
}

temp2$avg_value2 <- 90 - temp2$avg_value
temp2$avg_value2[which(temp2$avg_value2 < 0)] <- temp2$avg_value2[which(temp2$avg_value2 < 0)] + 360
favstats(temp2$avg_value2~temp2$rainy)
## dry: 281.70 --> WNW
## rainy: 6.11 --> N

#
#
### Time at locations ####
clust_time_part_median_dry <- c((favstats(Choma5_rain_clusts_shorter2_no_home_dry$clust_trips_time_new~Choma5_rain_clusts_shorter2_no_home_dry$partid)$median)*24,
                                   (favstats(Mutasa5_rain_clusts_shorter2_no_home_dry$clust_trips_time_new~Mutasa5_rain_clusts_shorter2_no_home_dry$partid)$median)*24,
                                   (favstats(Nchelenge5_clusts_shorter2_no_home_dry$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_dry$partid)$median)*24)
clust_time_part_median_rainy <- c((favstats(Choma5_rain_clusts_shorter2_no_home_rainy$clust_trips_time_new~Choma5_rain_clusts_shorter2_no_home_rainy$partid)$median)*24,
                                 (favstats(Mutasa5_rain_clusts_shorter2_no_home_rainy$clust_trips_time_new~Mutasa5_rain_clusts_shorter2_no_home_rainy$partid)$median)*24,
                                 (favstats(Nchelenge5_clusts_shorter2_no_home_rainy$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_rainy$partid)$median)*24)
favstats(clust_time_part_median_dry)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4605556 2.099028 3.232778 4.821944 195.3874 6.545787 19.76308 133       0
favstats(clust_time_part_median_rainy)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.8686111 1.850833 3.085417 4.335556 98.0725 5.712681 12.09442 139       0
clust_time_part_median_dry2 <- clust_time_part_median_dry[which(clust_time_part_median_dry < 35)]
clust_time_part_median_rainy2 <- clust_time_part_median_rainy[which(clust_time_part_median_rainy < 35)]
favstats(clust_time_part_median_dry2)
#       min       Q1   median       Q3      max     mean      sd  n missing
# 0.4605556 2.072465 3.213472 4.613715 19.32806 3.860087 2.759338 130       0
favstats(clust_time_part_median_rainy2)
#       min        Q1   median       Q3      max     mean       sd  n missing
# 0.8686111 1.820278 3.0175 4.193403 22.93264 3.856413 3.382445 135       0

clust_time_part_gender <- as.data.frame(cbind(c(clust_time_part_median_dry2, clust_time_part_median_rainy2),
                                              c(rep("dry",length(clust_time_part_median_dry2)),
                                                rep("rainy",length(clust_time_part_median_rainy2)))))
colnames(clust_time_part_gender) <- c("values","type")
clust_time_part_gender$values <- as.numeric(as.character(clust_time_part_gender$values))
clust_time_part_gender$type <- factor(clust_time_part_gender$type)

ggplot(data=clust_time_part_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,24)) +
  scale_x_continuous(breaks=c(seq(0,30,2)))


shapiro.test(clust_time_part_gender$values) ## not normal
kruskal.test(clust_time_part_gender$values, clust_time_part_gender$type) # p = 0.3

my_comp <- list(c("rainy","dry"))
ggplot(data=clust_time_part_gender, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  stat_compare_means(comparisons = my_comp,method="wilcox.test", label="p.signif", size=5.5, color="red", hide.ns=TRUE) +
  stat_compare_means(label.y = 25) 

#
#
#### Biting time metrics #####
### number locations ####
loc_counts_dry <- c(favstats(Choma5_clusts_shorter2_no_home_biting_places2_dry$clust~Choma5_clusts_shorter2_no_home_biting_places2_dry$partid)$n,
                       favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_dry$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_dry$partid)$n,
                       favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_dry$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_dry$partid)$n)
loc_counts_rainy <- c(favstats(Choma5_clusts_shorter2_no_home_biting_places2_rainy$clust~Choma5_clusts_shorter2_no_home_biting_places2_rainy$partid)$n,
                     favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_rainy$clust~Mutasa5_clusts_shorter2_no_home_biting_places2_rainy$partid)$n,
                     favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy$clust~Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy$partid)$n)
only_home_dry1 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_dry$partid))
only_home_rainy1 <- nlevels(as.factor(Choma5_clusts_shorter2$partid[which(Choma5_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_biting_places2_rainy$partid))
only_home_dry2 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_dry$partid))
only_home_rainy2 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid[which(Mutasa5_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_biting_places2_rainy$partid))
only_home_dry3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$rainy == 0)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_dry$partid))
only_home_rainy3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid[which(Nchelenge5_clusts_shorter2$rainy == 1)])) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy$partid))

loc_counts_dry2 <- c(loc_counts_dry,rep(0, only_home_dry1+only_home_dry2+only_home_dry3))
loc_counts_rainy2 <- c(loc_counts_rainy,rep(0, only_home_rainy1 + only_home_rainy2 + only_home_rainy3))

locs_count2 <- as.data.frame(cbind(c(loc_counts_dry2, loc_counts_rainy2),
                                   c(rep("dry",length(loc_counts_dry2)),
                                     rep("rainy",length(loc_counts_rainy2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type)

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) + 
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(-0.5,13)) +
  scale_x_continuous(breaks=c(seq(0,25,2)))
##
### time per location ####
clust_time_part_median_dry <- c((favstats(Choma5_clusts_shorter2_no_home_biting_places2_dry$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_dry$partid)$median)*24,
                                   (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_dry$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_dry$partid)$median)*24,
                                   (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_dry$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_dry$partid)$median)*24)
clust_time_part_median_rainy <- c((favstats(Choma5_clusts_shorter2_no_home_biting_places2_rainy$clust_trips_time_new_biting_time~Choma5_clusts_shorter2_no_home_biting_places2_rainy$partid)$median)*24,
                                 (favstats(Mutasa5_clusts_shorter2_no_home_biting_places2_rainy$clust_trips_time_new_biting_time~Mutasa5_clusts_shorter2_no_home_biting_places2_rainy$partid)$median)*24,
                                 (favstats(Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy$clust_trips_time_new_biting_time~Nchelenge5_clusts_shorter2_no_home_biting_places2_rainy$partid)$median)*24)
clust_time_part_median_rainy2 <- clust_time_part_median_rainy[which(clust_time_part_median_rainy < 500)]

clust_time_part_gender <- as.data.frame(cbind(c(clust_time_part_median_dry, clust_time_part_median_rainy2),
                                              c(rep("dry",length(clust_time_part_median_dry)),
                                                rep("rainy",length(clust_time_part_median_rainy2)))))
colnames(clust_time_part_gender) <- c("values","type")
clust_time_part_gender$values <- as.numeric(as.character(clust_time_part_gender$values))
clust_time_part_gender$type <- factor(clust_time_part_gender$type)

ggplot(data=clust_time_part_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=4) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,300)) +
  scale_x_continuous(breaks=c(seq(0,300,50)))
##
### number trips ####
### PER PART ###
trips_part_median_dry <- c(favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_dry$partid)$median,
                              favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_dry$partid)$median,
                              favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_dry$partid)$median)
trips_part_median_rainy <- c(favstats(Choma5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2~Choma5_clust_trips_shorter2_no_home_biting_places3_rainy$partid)$median,
                            favstats(Mutasa5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2~Mutasa5_clust_trips_shorter2_no_home_biting_places3_rainy$partid)$median,
                            favstats(Nchelenge5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2~Nchelenge5_clust_trips_shorter2_no_home_biting_places3_rainy$partid)$median)
trips_part_median_dry2 <- trips_part_median_dry[which(trips_part_median_dry < 20)]


trips_part2_gender <- as.data.frame(cbind(c(trips_part_median_dry2, trips_part_median_rainy),
                                          c(rep("dry",length(trips_part_median_dry2)),
                                            rep("rainy",length(trips_part_median_rainy)))))
colnames(trips_part2_gender) <- c("values","type")
trips_part2_gender$values <- as.numeric(as.character(trips_part2_gender$values))
trips_part2_gender$type <- factor(trips_part2_gender$type)

ggplot(data=trips_part2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,14)) + 
  scale_x_continuous(breaks=c(seq(0,14,2)))



dry_all <- c(Choma5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2, 
                Mutasa5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2,
                Nchelenge5_clust_trips_shorter2_no_home_biting_places3_dry$trip_count_new2)
rainy_all <- c(Choma5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2, 
              Mutasa5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2,
              Nchelenge5_clust_trips_shorter2_no_home_biting_places3_rainy$trip_count_new2)

trips_all2 <- as.data.frame(cbind(c(dry_all, rainy_all),
                                  c(rep("dry",length(dry_all)),
                                    rep("rainy",length(rainy_all)))))
colnames(trips_all2) <- c("values","type")
trips_all2$values <- as.numeric(as.character(trips_all2$values))
trips_all2$type <- factor(trips_all2$type)

table(trips_all2$type, as.factor(trips_all2$values))
#
##### percent biting time spent at home vs elsewhere #####
perc_dry <- c(Choma5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time, 
                 Mutasa5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time, 
                 Nchelenge5_clusts_shorter2_biting_places3b_dry$percent_home_biting_time)
perc_rainy <- c(Choma5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time, 
               Mutasa5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time, 
               Nchelenge5_clusts_shorter2_biting_places3b_rainy$percent_home_biting_time)

perc_time_home2_gender <- as.data.frame(cbind(c(perc_dry, perc_rainy),
                                              c(rep("dry",length(perc_dry)),
                                                rep("rainy",length(perc_rainy)))))
colnames(perc_time_home2_gender) <- c("values","type")
perc_time_home2_gender$values <- as.numeric(as.character(perc_time_home2_gender$values))
perc_time_home2_gender$type <- factor(perc_time_home2_gender$type)

ggplot(data=perc_time_home2_gender, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=1) +
  scale_fill_manual(values = c("sienna", "skyblue")) +
  scale_color_manual(values = c("sienna", "skyblue")) +
  coord_cartesian(xlim=c(0,100))
##
#############################################
#############################################
##### compare between time of day without locations as covariate #####
### number of locations (without home) ####
# -Number of locations where at least 1/4 time spent there during “time of day”

loc_counts_morning <- c(favstats(Choma5_clusts_shorter2_no_home_morning$clust~Choma5_clusts_shorter2_no_home_morning$partid)$n,
                        favstats(Mutasa5_clusts_shorter2_no_home_morning$clust~Mutasa5_clusts_shorter2_no_home_morning$partid)$n,
                        favstats(Nchelenge5_clusts_shorter2_no_home_morning$clust~Nchelenge5_clusts_shorter2_no_home_morning$partid)$n)
loc_counts_midday <- c(favstats(Choma5_clusts_shorter2_no_home_midday$clust~Choma5_clusts_shorter2_no_home_midday$partid)$n,
                        favstats(Mutasa5_clusts_shorter2_no_home_midday$clust~Mutasa5_clusts_shorter2_no_home_midday$partid)$n,
                        favstats(Nchelenge5_clusts_shorter2_no_home_midday$clust~Nchelenge5_clusts_shorter2_no_home_midday$partid)$n)
loc_counts_evening <- c(favstats(Choma5_clusts_shorter2_no_home_evening$clust~Choma5_clusts_shorter2_no_home_evening$partid)$n,
                        favstats(Mutasa5_clusts_shorter2_no_home_evening$clust~Mutasa5_clusts_shorter2_no_home_evening$partid)$n,
                        favstats(Nchelenge5_clusts_shorter2_no_home_evening$clust~Nchelenge5_clusts_shorter2_no_home_evening$partid)$n)
loc_counts_night <- c(favstats(Choma5_clusts_shorter2_no_home_night$clust~Choma5_clusts_shorter2_no_home_night$partid)$n,
                        favstats(Mutasa5_clusts_shorter2_no_home_night$clust~Mutasa5_clusts_shorter2_no_home_night$partid)$n,
                        favstats(Nchelenge5_clusts_shorter2_no_home_night$clust~Nchelenge5_clusts_shorter2_no_home_night$partid)$n)

only_home_morning1 <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_morning$partid))
only_home_midday1 <-  nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_midday$partid))
only_home_evening1 <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_evening$partid))
only_home_night1 <- nlevels(as.factor(Choma5_clusts_shorter2$partid)) - nlevels(as.factor(Choma5_clusts_shorter2_no_home_night$partid))

only_home_morning2 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_morning$partid))
only_home_midday2 <-  nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_midday$partid))
only_home_evening2 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_evening$partid))
only_home_night2 <- nlevels(as.factor(Mutasa5_clusts_shorter2$partid)) - nlevels(as.factor(Mutasa5_clusts_shorter2_no_home_night$partid))

only_home_morning3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_morning$partid))
only_home_midday3 <-  nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_midday$partid))
only_home_evening3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_evening$partid))
only_home_night3 <- nlevels(as.factor(Nchelenge5_clusts_shorter2$partid)) - nlevels(as.factor(Nchelenge5_clusts_shorter2_no_home_night$partid))


loc_counts_morning2 <- c(loc_counts_morning,rep(0, only_home_morning1+only_home_morning2+only_home_morning3))
loc_counts_midday2 <- c(loc_counts_midday,rep(0, only_home_midday1+only_home_midday2+only_home_midday3))
loc_counts_evening2 <- c(loc_counts_evening,rep(0, only_home_evening1+only_home_evening2+only_home_evening3))
loc_counts_night2 <- c(loc_counts_night,rep(0, only_home_night1+only_home_night2+only_home_night3))

locs_count2 <- as.data.frame(cbind(c(loc_counts_morning2, loc_counts_midday2, loc_counts_evening2, loc_counts_night2),
                                   c(rep("morning",length(loc_counts_morning2)),
                                     rep("midday",length(loc_counts_midday2)),
                                     rep("evening",length(loc_counts_evening2)),
                                     rep("night",length(loc_counts_night2)))))
colnames(locs_count2) <- c("values","type")
locs_count2$values <- as.numeric(as.character(locs_count2$values))
locs_count2$type <- factor(locs_count2$type, levels=c("morning","midday","evening","night"))
timecols <- c("#f8f871","#35c7fc", "#fC6a35", "#6a35fc")

ggplot(data=locs_count2, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 2, alpha=0.8) +
  geom_density(aes(color=type), adjust=3) + 
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-0.5,25)) +
  scale_x_continuous(breaks=c(seq(0,26,2)))


shapiro.test(locs_count2$values) ## not normal
kruskal.test(locs_count2$values, locs_count2$type) # p < 2e-16
pairwise.wilcox.test(locs_count2$values, locs_count2$type, p.adjust.method = "holm") ## 
#         morning midday evening
# midday  <2e-16     -      -      
# evening    0.5  <2e-16    -      
# night   <2e-16  <2e-16 <2e-16 


my_comp <- list(c("morning", "midday"), c("morning","evening"), c("morning","night"), c("midday","evening"), c("midday","night"), c("evening","night"))
ggplot(data=locs_count2, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_x_discrete("Time of day", labels=c("Morning", "Midday", "Evening", "Night")) +
  scale_y_continuous("Number of lcoations visited", breaks=c(seq(from=0,to=20,by=5))) +
  scale_fill_manual(values = timecols) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.8,step.increase=0.07 ) +
  guides(fill="none") +
  theme(axis.title=element_text(size=18),
        axis.text = element_text(size=16))

#
### distance of locations (without home) ####
# -Median distance of locations where at least 1/4 time spent there during “time of day”
medians_km_morning <- c((favstats(Choma5_clusts_shorter2_no_home_morning$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_morning$partid)$median),
                        (favstats(Mutasa5_clusts_shorter2_no_home_morning$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_morning$partid)$median),
                        (favstats(Nchelenge5_clusts_shorter2_no_home_morning$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_morning$partid)$median))
medians_km_midday <- c((favstats(Choma5_clusts_shorter2_no_home_midday$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_midday$partid)$median),
                        (favstats(Mutasa5_clusts_shorter2_no_home_midday$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_midday$partid)$median),
                        (favstats(Nchelenge5_clusts_shorter2_no_home_midday$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_midday$partid)$median))
medians_km_evening <- c((favstats(Choma5_clusts_shorter2_no_home_evening$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_evening$partid)$median),
                        (favstats(Mutasa5_clusts_shorter2_no_home_evening$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_evening$partid)$median),
                        (favstats(Nchelenge5_clusts_shorter2_no_home_evening$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_evening$partid)$median))
medians_km_night <- c((favstats(Choma5_clusts_shorter2_no_home_night$hhdist_km_haversine~Choma5_clusts_shorter2_no_home_night$partid)$median),
                        (favstats(Mutasa5_clusts_shorter2_no_home_night$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home_night$partid)$median),
                        (favstats(Nchelenge5_clusts_shorter2_no_home_night$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home_night$partid)$median))
medians_km_morning_short <- medians_km_morning[which(medians_km_morning < 240)]
medians_km_midday_short <- medians_km_midday[which(medians_km_midday < 240)]
medians_km_evening_short <- medians_km_evening[which(medians_km_evening < 240)]
medians_km_night_short <- medians_km_night[which(medians_km_night < 240)]

dists_km2_time <- as.data.frame(cbind(c(medians_km_morning_short, medians_km_midday_short, medians_km_evening_short, medians_km_night_short),
                                      c(rep("morning",length(medians_km_morning_short)),
                                        rep("midday",length(medians_km_midday_short)),
                                        rep("evening",length(medians_km_evening_short)),
                                        rep("night",length(medians_km_night_short)))))
colnames(dists_km2_time) <- c("values","type")
dists_km2_time$values <- as.numeric(as.character(dists_km2_time$values))
dists_km2_time$type <- factor(dists_km2_time$type, levels=c("morning","midday","evening","night"))

ggplot(data=dists_km2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 25, alpha=0.8) +
  geom_density(aes(color=type), adjust=15) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-4,245))


shapiro.test(dists_km2_time$values) ## not normal
kruskal.test(dists_km2_time$values, dists_km2_time$type) # p = 9e-06
pairwise.wilcox.test(dists_km2_time$values, dists_km2_time$type, p.adjust.method = "holm") ## 
#         morning midday evening
# midday  0.9     -      -      
# evening    0.9  0.9     -      
# night   1e-04  3e-06 1e-04


my_comp <- list(c("morning", "midday"), c("morning","evening"), c("morning","night"), c("midday","evening"), c("midday","night"), c("evening","night"))
ggplot(data=dists_km2_time, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = timecols) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.7) +
  stat_compare_means(label.y = 385) 

#
###
### number of trips (without home) ####
# -Median number of trips where at least 1/4 of the trip time was during “time of day”
trips_part_median_morning <- c(favstats(Choma_number_trips_morning$trips ~ Choma_number_trips_morning$partid)$median,
                               favstats(Mutasa_number_trips_morning$trips ~ Mutasa_number_trips_morning$partid)$median,
                               favstats(Nchelenge_number_trips_morning$trips ~ Nchelenge_number_trips_morning$partid)$median)
trips_part_median_midday <- c(favstats(Choma_number_trips_midday$trips ~ Choma_number_trips_midday$partid)$median,
                               favstats(Mutasa_number_trips_midday$trips ~ Mutasa_number_trips_midday$partid)$median,
                               favstats(Nchelenge_number_trips_midday$trips ~ Nchelenge_number_trips_midday$partid)$median)
trips_part_median_evening <- c(favstats(Choma_number_trips_evening$trips ~ Choma_number_trips_evening$partid)$median,
                               favstats(Mutasa_number_trips_evening$trips ~ Mutasa_number_trips_evening$partid)$median,
                               favstats(Nchelenge_number_trips_evening$trips ~ Nchelenge_number_trips_evening$partid)$median)
trips_part_median_night <- c(favstats(Choma_number_trips_night$trips ~ Choma_number_trips_night$partid)$median,
                               favstats(Mutasa_number_trips_night$trips ~ Mutasa_number_trips_night$partid)$median,
                               favstats(Nchelenge_number_trips_night$trips ~ Nchelenge_number_trips_night$partid)$median)

trips_part2_time <- as.data.frame(cbind(c(trips_part_median_morning, trips_part_median_midday, trips_part_median_evening, trips_part_median_night),
                                        c(rep("morning",length(trips_part_median_morning)),
                                          rep("midday",length(trips_part_median_midday)),
                                          rep("evening",length(trips_part_median_evening)),
                                          rep("night",length(trips_part_median_night)))))
colnames(trips_part2_time) <- c("values","type")
trips_part2_time$values <- as.numeric(as.character(trips_part2_time$values))
trips_part2_time$type <- factor(trips_part2_time$type, levels=c("morning","midday","evening","night"))

ggplot(data=trips_part2_time[which(trips_part2_time$values < 13),], aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 1, alpha=0.8) +
  geom_density(aes(color=type), adjust=2.5) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(0.5,16.25)) +
  scale_x_continuous(breaks=c(seq(0,20,2)))


shapiro.test(trips_part2_time$values) ## not normal
kruskal.test(trips_part2_time$values, trips_part2_time$type) # p = 0.05
pairwise.wilcox.test(trips_part2_time$values, trips_part2_time$type, p.adjust.method = "holm") ## 
#         morning midday evening
# midday  1     -      -      
# evening   1  1     -      
# night   0.1  0.05 0.12


my_comp <- list(c("morning", "midday"), c("morning","evening"), c("morning","night"), c("midday","evening"), c("midday","night"), c("evening","night"))
ggplot(data=trips_part2_time, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = timecols) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.7) +
  stat_compare_means(label.y = 27) 

#
### time per trip ####
# -Median time per trip for trips where at least 1/4 trip time was during “time of day”
trip_time_part_median_morning <- c(favstats(Choma5_clust_trips_shorter2_no_home_morning$trip_time_new2 ~ Choma5_clust_trips_shorter2_no_home_morning$partid)$median,
                                   favstats(Mutasa5_clust_trips_shorter2_no_home_morning$trip_time_new2 ~ Mutasa5_clust_trips_shorter2_no_home_morning$partid)$median,
                                   favstats(Nchelenge5_clust_trips_shorter2_no_home_morning$trip_time_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_morning$partid)$median)
trip_time_part_median_midday <- c(favstats(Choma5_clust_trips_shorter2_no_home_midday$trip_time_new2 ~ Choma5_clust_trips_shorter2_no_home_midday$partid)$median,
                                   favstats(Mutasa5_clust_trips_shorter2_no_home_midday$trip_time_new2 ~ Mutasa5_clust_trips_shorter2_no_home_midday$partid)$median,
                                   favstats(Nchelenge5_clust_trips_shorter2_no_home_midday$trip_time_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_midday$partid)$median)
trip_time_part_median_evening <- c(favstats(Choma5_clust_trips_shorter2_no_home_evening$trip_time_new2 ~ Choma5_clust_trips_shorter2_no_home_evening$partid)$median,
                                   favstats(Mutasa5_clust_trips_shorter2_no_home_evening$trip_time_new2 ~ Mutasa5_clust_trips_shorter2_no_home_evening$partid)$median,
                                   favstats(Nchelenge5_clust_trips_shorter2_no_home_evening$trip_time_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_evening$partid)$median)
trip_time_part_median_night <- c(favstats(Choma5_clust_trips_shorter2_no_home_night$trip_time_new2 ~ Choma5_clust_trips_shorter2_no_home_night$partid)$median,
                                   favstats(Mutasa5_clust_trips_shorter2_no_home_night$trip_time_new2 ~ Mutasa5_clust_trips_shorter2_no_home_night$partid)$median,
                                   favstats(Nchelenge5_clust_trips_shorter2_no_home_night$trip_time_new2 ~ Nchelenge5_clust_trips_shorter2_no_home_night$partid)$median)


trip_time_part_median_morning_short <- trip_time_part_median_morning[which(trip_time_part_median_morning*24 < 150)]
trip_time_part_median_midday_short <- trip_time_part_median_midday[which(trip_time_part_median_midday*24 < 150)]
trip_time_part_median_evening_short <- trip_time_part_median_evening[which(trip_time_part_median_evening*24 < 150)] 
trip_time_part_median_night_short <- trip_time_part_median_night[which(trip_time_part_median_night*24 < 150)] 

trip_times_part2_time <- as.data.frame(cbind(c(trip_time_part_median_morning_short*24, trip_time_part_median_midday_short*24, trip_time_part_median_evening_short*24, trip_time_part_median_night_short*24),
                                             c(rep("morning",length(trip_time_part_median_morning_short)),
                                               rep("midday",length(trip_time_part_median_midday_short)),
                                               rep("evening",length(trip_time_part_median_evening_short)),
                                               rep("night",length(trip_time_part_median_night_short)))))
colnames(trip_times_part2_time) <- c("values","type")
trip_times_part2_time$values <- as.numeric(as.character(trip_times_part2_time$values))
trip_times_part2_time$type <- factor(trip_times_part2_time$type, levels=c("morning","midday","evening","night"))


ggplot(data=trip_times_part2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 10, alpha=0.8) +
  geom_density(aes(color=type), adjust=9) +
  scale_fill_manual("Time of day", values = timecols,  labels=c("Morning", "Midday", "Evening", "Night")) +
  scale_color_manual("Time of day", values = timecols, labels=c("Morning", "Midday", "Evening", "Night")) +
  coord_cartesian(xlim=c(-2,130)) + 
  scale_x_continuous("Median time per trip (hr)", breaks=c(seq(0,135,20))) +
  theme_classic() +
  theme(axis.title=element_text(size=14),
        axis.text = element_text(size=12),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.position=c(0.88,0.3))


shapiro.test(trip_times_part2_time$values) ## not normal
kruskal.test(trip_times_part2_time$values, trip_times_part2_time$type) # p < 2e-16
pairwise.wilcox.test(trip_times_part2_time$values, trip_times_part2_time$type, p.adjust.method = "holm") ## 
#         morning midday evening
#   midday  0.003   -      -      
#   evening 1e-05   0.011  -      
#   night  <2e-16  <2e-16 <2e-16 


my_comp <- list(c("morning", "midday"), c("morning","evening"), c("morning","night"), c("midday","evening"), c("midday","night"), c("evening","night"))
ggplot(data=trip_times_part2_time, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = timecols) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.7) +
  stat_compare_means(label.y = 230) 

#
### Percent time at home vs elsewhere ####
# -Percent time at home (vs elsewhere) during “time of day”
perc_home_morning <- c((Choma5_clusts_shorter2_home$clust_trips_time_new_morning/Choma5_clusts_shorter2_home$part_trips_time_new_morning)*100,
                       (Mutasa5_clusts_shorter2_home$clust_trips_time_new_morning/Mutasa5_clusts_shorter2_home$part_trips_time_new_morning)*100,
                       (Nchelenge5_clusts_shorter2_home$clust_trips_time_new_morning/Nchelenge5_clusts_shorter2_home$part_trips_time_new_morning)*100)
perc_home_midday <- c((Choma5_clusts_shorter2_home$clust_trips_time_new_midday/Choma5_clusts_shorter2_home$part_trips_time_new_midday)*100,
                       (Mutasa5_clusts_shorter2_home$clust_trips_time_new_midday/Mutasa5_clusts_shorter2_home$part_trips_time_new_midday)*100,
                       (Nchelenge5_clusts_shorter2_home$clust_trips_time_new_midday/Nchelenge5_clusts_shorter2_home$part_trips_time_new_midday)*100)
perc_home_evening <- c((Choma5_clusts_shorter2_home$clust_trips_time_new_evening/Choma5_clusts_shorter2_home$part_trips_time_new_evening)*100,
                       (Mutasa5_clusts_shorter2_home$clust_trips_time_new_evening/Mutasa5_clusts_shorter2_home$part_trips_time_new_evening)*100,
                       (Nchelenge5_clusts_shorter2_home$clust_trips_time_new_evening/Nchelenge5_clusts_shorter2_home$part_trips_time_new_evening)*100)
perc_home_night <- c((Choma5_clusts_shorter2_home$clust_trips_time_new_night/Choma5_clusts_shorter2_home$part_trips_time_new_night)*100,
                       (Mutasa5_clusts_shorter2_home$clust_trips_time_new_night/Mutasa5_clusts_shorter2_home$part_trips_time_new_night)*100,
                       (Nchelenge5_clusts_shorter2_home$clust_trips_time_new_night/Nchelenge5_clusts_shorter2_home$part_trips_time_new_night)*100)

perc_home_morning[which(is.na(perc_home_morning))] <- 0
perc_home_midday[which(is.na(perc_home_midday))] <- 0
perc_home_evening[which(is.na(perc_home_evening))] <- 0
perc_home_night[which(is.na(perc_home_night))] <- 0

perc_home2_time <- as.data.frame(cbind(c(perc_home_morning, perc_home_midday, perc_home_evening, perc_home_night),
                                       c(rep("morning",length(perc_home_morning)),
                                         rep("midday",length(perc_home_midday)),
                                         rep("evening",length(perc_home_evening)),
                                         rep("night",length(perc_home_night)))))
colnames(perc_home2_time) <- c("values","type")
perc_home2_time$values <- as.numeric(as.character(perc_home2_time$values))
perc_home2_time$type <- factor(perc_home2_time$type, levels=c("morning","midday","evening","night"))

ggplot(data=perc_home2_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 20, alpha=0.8) +
  geom_density(aes(color=type), adjust=1.5) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-7,107)) + 
  scale_x_continuous(breaks=c(seq(0,100,20)))



shapiro.test(perc_home2_time$values) ## not normal
kruskal.test(perc_home2_time$values, perc_home2_time$type) # p < 2e-16
pairwise.wilcox.test(perc_home2_time$values, perc_home2_time$type, p.adjust.method = "holm") ## 
#         morning midday evening
#   midday  2e-08   -      -      
#   evening 0.1     3e-11  -      
#   night   4e-06   2e-15  7e-05  


my_comp <- list(c("morning", "midday"), c("morning","evening"), c("morning","night"), c("midday","evening"), c("midday","night"), c("evening","night"))
ggplot(data=perc_home2_time, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = timecols) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.7) +
  stat_compare_means(label.y = 170) 

##
### Time at locations ####
clust_time_part_median_morning <- c((favstats(Choma5_clusts_shorter2_no_home_morning$clust_trips_time_new~Choma5_clusts_shorter2_no_home_morning$partid)$median)*24,
                                    (favstats(Mutasa5_clusts_shorter2_no_home_morning$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_morning$partid)$median)*24,
                                    (favstats(Nchelenge5_clusts_shorter2_no_home_morning$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_morning$partid)$median)*24)
clust_time_part_median_midday <- c((favstats(Choma5_clusts_shorter2_no_home_midday$clust_trips_time_new~Choma5_clusts_shorter2_no_home_midday$partid)$median)*24,
                                    (favstats(Mutasa5_clusts_shorter2_no_home_midday$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_midday$partid)$median)*24,
                                    (favstats(Nchelenge5_clusts_shorter2_no_home_midday$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_midday$partid)$median)*24)
clust_time_part_median_evening <- c((favstats(Choma5_clusts_shorter2_no_home_evening$clust_trips_time_new~Choma5_clusts_shorter2_no_home_evening$partid)$median)*24,
                                    (favstats(Mutasa5_clusts_shorter2_no_home_evening$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_evening$partid)$median)*24,
                                    (favstats(Nchelenge5_clusts_shorter2_no_home_evening$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_evening$partid)$median)*24)
clust_time_part_median_night <- c((favstats(Choma5_clusts_shorter2_no_home_night$clust_trips_time_new~Choma5_clusts_shorter2_no_home_night$partid)$median)*24,
                                    (favstats(Mutasa5_clusts_shorter2_no_home_night$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home_night$partid)$median)*24,
                                    (favstats(Nchelenge5_clusts_shorter2_no_home_night$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home_night$partid)$median)*24)

clust_time_part_median_morning2 <- clust_time_part_median_morning[which(clust_time_part_median_morning < 400)]
clust_time_part_median_midday2 <- clust_time_part_median_midday[which(clust_time_part_median_midday < 400)]
clust_time_part_median_evening2 <- clust_time_part_median_evening[which(clust_time_part_median_evening < 400)]
clust_time_part_median_night2 <- clust_time_part_median_night[which(clust_time_part_median_night < 400)]

clust_time_part_time <- as.data.frame(cbind(c(clust_time_part_median_morning2, clust_time_part_median_midday2,clust_time_part_median_evening2,clust_time_part_median_night2),
                                            c(rep("morning",length(clust_time_part_median_morning2)),
                                              rep("midday",length(clust_time_part_median_midday2)),
                                              rep("evening",length(clust_time_part_median_evening2)),
                                              rep("night",length(clust_time_part_median_night2)))))
colnames(clust_time_part_time) <- c("values","type")
clust_time_part_time$values <- as.numeric(as.character(clust_time_part_time$values))
clust_time_part_time$type <- factor(clust_time_part_time$type, levels=c("morning","midday","evening","night"))

ggplot(data=clust_time_part_time, aes(x=values, y=..density..)) +
  geom_histogram(aes(fill=type), position = "dodge", binwidth = 50, alpha=0.8) +
  geom_density(aes(color=type), adjust=25) +
  scale_fill_manual(values = timecols) +
  scale_color_manual(values = timecols) +
  coord_cartesian(xlim=c(-10,400)) +
  scale_x_continuous(breaks=c(seq(0,900,50)))


shapiro.test(clust_time_part_time$values) ## not normal
kruskal.test(clust_time_part_time$values, clust_time_part_time$type) # p < 2e-16
pairwise.wilcox.test(clust_time_part_time$values, clust_time_part_time$type, p.adjust.method = "holm") ## 
#         morning midday evening
#   midday  1e-05   -      -      
#   evening 0.001   0.764  -      
#   night   <2e-16  <2e-16 <2e-16 


my_comp <- list(c("morning", "midday"), c("morning","evening"), c("morning","night"), c("midday","evening"), c("midday","night"), c("evening","night"))
ggplot(data=clust_time_part_time, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = timecols) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.7) +
  stat_compare_means(label.y = 650) 

#
#
#
### Average direction ####
######### all ###
Cho_ind <- which(colnames(Choma5_clusts_shorter2_no_home_night) %in% c("partid", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))
Mut_ind <- which(colnames(Mutasa5_clusts_shorter2_no_home_night) %in% c("partid", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))
Nch_ind <- which(colnames(Nchelenge5_clusts_shorter2_no_home_night) %in% c("partid", "avg_dir_from_home_bearing","avg_dir_from_home_mag", "avg_dir_from_home_values"))

temp_morning <- rbind(Choma5_clusts_shorter2_no_home_morning[!duplicated(Choma5_clusts_shorter2_no_home_morning$partid),c(Cho_ind)],
                  Mutasa5_clusts_shorter2_no_home_morning[!duplicated(Mutasa5_clusts_shorter2_no_home_morning$partid),c(Mut_ind)],
                  Nchelenge5_clusts_shorter2_no_home_morning[!duplicated(Nchelenge5_clusts_shorter2_no_home_morning$partid),c(Nch_ind)])
temp_midday <- rbind(Choma5_clusts_shorter2_no_home_midday[!duplicated(Choma5_clusts_shorter2_no_home_midday$partid),c(Cho_ind)],
                  Mutasa5_clusts_shorter2_no_home_midday[!duplicated(Mutasa5_clusts_shorter2_no_home_midday$partid),c(Mut_ind)],
                  Nchelenge5_clusts_shorter2_no_home_midday[!duplicated(Nchelenge5_clusts_shorter2_no_home_midday$partid),c(Nch_ind)])
temp_evening <- rbind(Choma5_clusts_shorter2_no_home_evening[!duplicated(Choma5_clusts_shorter2_no_home_evening$partid),c(Cho_ind)],
                  Mutasa5_clusts_shorter2_no_home_evening[!duplicated(Mutasa5_clusts_shorter2_no_home_evening$partid),c(Mut_ind)],
                  Nchelenge5_clusts_shorter2_no_home_evening[!duplicated(Nchelenge5_clusts_shorter2_no_home_evening$partid),c(Nch_ind)])
temp_night <- rbind(Choma5_clusts_shorter2_no_home_night[!duplicated(Choma5_clusts_shorter2_no_home_night$partid),c(Cho_ind)],
                  Mutasa5_clusts_shorter2_no_home_night[!duplicated(Mutasa5_clusts_shorter2_no_home_night$partid),c(Mut_ind)],
                  Nchelenge5_clusts_shorter2_no_home_night[!duplicated(Nchelenge5_clusts_shorter2_no_home_night$partid),c(Nch_ind)])


temp_dir <- as.data.frame(cbind(c(temp_morning$avg_dir_from_home_bearing, temp_midday$avg_dir_from_home_bearing,temp_evening$avg_dir_from_home_bearing, temp_night$avg_dir_from_home_bearing),
                                c(temp_morning$avg_dir_from_home_mag, temp_midday$avg_dir_from_home_mag,temp_evening$avg_dir_from_home_mag, temp_night$avg_dir_from_home_mag),
                                c(temp_morning$avg_dir_from_home_values, temp_midday$avg_dir_from_home_values,temp_evening$avg_dir_from_home_values, temp_night$avg_dir_from_home_values),
                                c(rep("morning",nrow(temp_morning)),
                                  rep("midday",nrow(temp_midday)),
                                  rep("evening",nrow(temp_evening)),
                                  rep("night",nrow(temp_night)))))

colnames(temp_dir) <- c("bearing","mag", "values", "type")
temp_dir$bearing <- as.numeric(as.character(temp_dir$bearing))
temp_dir$mag <- as.numeric(as.character(temp_dir$mag))
temp_dir$values <- as.numeric(as.character(temp_dir$values))
temp_dir$type <- factor(temp_dir$type, levels=c("morning","midday","evening","night"))


ggplot(data=temp_dir, aes(x=bearing, y=mag)) +
  geom_col(aes(fill=type), position=position_dodge2(preserve = "single"), width=1) +
  coord_polar(theta="x",start=0) +
  scale_x_continuous("Average direction of locations", limits=c(0,360),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                               "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  scale_y_continuous("Average distance of locations (km)") +
  theme_bw() +
  scale_fill_manual("Time of Day", values=timecols, labels=c("morning", "midday", "evening", "night"))


temp_dir$values[which(temp_dir$values > 348.75)] <- -(360-temp_dir$values[which(temp_dir$values > 348.75)])
ggplot(data=temp_dir, aes(x=values, y=..count.., fill=type)) +
  geom_histogram(position = "dodge", breaks=c(seq(-11.25,348.75,22.5)), alpha=0.8) +
  coord_polar(theta="x",start=(-pi/16)) +
  scale_x_continuous("Average direction of locations", limits=c(-11.25,348.75),breaks=c(seq(0,338.5,22.5)), labels = c("N", "NNE", "NE", "ENE","E", "ESE", "SE", "SSE",
                                                                                                                       "S", "SSW", "SW", "WSW","W", "WNW", "NW", "NNW")) +
  theme_bw() +
  theme(legend.position="none") +
  facet_wrap(~type)+
  scale_fill_manual(values=timecols) + 
  scale_y_continuous(breaks=c(seq(0,25,5)))




shapiro.test(temp_dir$values) ## not normal
kruskal.test(temp_dir$values, temp_dir$type) # p = 0.9
pairwise.wilcox.test(temp_dir$values, temp_dir$type, p.adjust.method = "holm") ## 
#         morning midday evening
#   midday  1   -      -      
#   evening 1   1  -      
#   night   1  1 1 


my_comp <- list(c("morning", "midday"), c("morning","evening"), c("morning","night"), c("midday","evening"), c("midday","night"), c("evening","night"))
ggplot(data=temp_dir, aes(x=type, y=values)) +
  geom_boxplot(aes(fill=type)) +  
  theme_classic() +
  scale_fill_manual(values = timecols) +
  stat_compare_means(comparisons = my_comp, label="p.signif", size=6, color="red", hide.ns=TRUE, vjust=0.7) +
  stat_compare_means(label.y = 600) 








## average direction calc
temp_dir$values2 <- 90 - temp_dir$values
temp_dir$values2[which(temp_dir$values2 < 0)] <- temp_dir$values2[which(temp_dir$values2 < 0)] + 360
temp_dir$values3 <- deg2rad(temp_dir$values2)

temp_dir$avg_value <- as.numeric(NA)
for(ii in 1:nlevels(temp_dir$type)){
  level <- levels(temp_dir$type)[ii]
  temp <- temp_dir[which(temp_dir$type == level),]
  sum_cos_dist <- sum(cos(temp$values3))/nrow(temp)
  sum_sin_dist <- sum(sin(temp$values3))/nrow(temp)
  avg_vec_dist <- sqrt((sum_cos_dist^2) + (sum_sin_dist^2))
  avg_deg <- rad2deg(atan2(sum_sin_dist,sum_cos_dist))
  if(avg_deg < 0){
    avg_deg <- 360+avg_deg
  }
  temp_dir$avg_value[which(temp_dir$type == level)] <- avg_deg
}

temp_dir$avg_value2 <- 90 - temp_dir$avg_value
temp_dir$avg_value2[which(temp_dir$avg_value2 < 0)] <- temp_dir$avg_value2[which(temp_dir$avg_value2 < 0)] + 360
favstats(temp_dir$avg_value2~temp_dir$type)
## morn: 263.00 --> W
## midd: 327.05 --> NNW
## even: 92.13 --> E
## nigh: 257.90 --> WSW
#
#
######
#############################################
#############################################
########
########
########
########
#############################################
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,9)]

loc_count <- as.data.frame(cbind(c(loc_counts_Choma$`Choma5_clusts_shorter2_no_home$partid`, 
                                       loc_counts_Mutasa$`Mutasa5_clusts_shorter2_no_home$partid`, 
                                       loc_counts_Nchelenge$`Nchelenge5_clusts_shorter2_no_home$partid`),
                                     c(loc_counts_Choma$n, 
                                       loc_counts_Mutasa$n, 
                                       loc_counts_Nchelenge$n),
                                     c(rep("Choma",length(loc_counts_Choma$n)),
                                       rep("Mutasa",length(loc_counts_Mutasa$n)),
                                       rep("Nchelenge",length(loc_counts_Nchelenge$n)))))
colnames(loc_count) <- c("partid","num_locs","loc")
loc_count$num_locs <- as.numeric(as.character(loc_count$num_locs))
loc_count$partid <- factor(loc_count$partid)
loc_count$loc <- factor(loc_count$loc)

loc_count$age <- as.numeric(NA)
loc_count$male <- factor(NA, levels=c(0,1))
loc_count$season <- factor(NA, levels=c("dry","rainy"))
for(ii in 1:length(loc_count$partid)){
  temp_part <- loc_count$partid[ii]
  a <- c(Choma5_clusts_shorter2_no_home$age_all[which(Choma5_clusts_shorter2_no_home$partid == temp_part)][1],
         Mutasa5_clusts_shorter2_no_home$age_all[which(Mutasa5_clusts_shorter2_no_home$partid == temp_part)][1],
         Nchelenge5_clusts_shorter2_no_home$age_all[which(Nchelenge5_clusts_shorter2_no_home$partid == temp_part)][1])
  a2 <- a[which(!(is.na(a)))]
  if(length(a2) > 0){
    loc_count$age[ii] <- a2
  }
  b <- c(as.character(Choma5_clusts_shorter2_no_home$male[which(Choma5_clusts_shorter2_no_home$partid == temp_part)][1]),
         as.character(Mutasa5_clusts_shorter2_no_home$male[which(Mutasa5_clusts_shorter2_no_home$partid == temp_part)][1]),
         as.character(Nchelenge5_clusts_shorter2_no_home$male[which(Nchelenge5_clusts_shorter2_no_home$partid == temp_part)][1]))
  b2 <- b[which(!(is.na(b)))]
  if(length(b2) > 0){
    loc_count$male[ii] <- b2
  }
  if(temp_part %in% c(unique(Choma5_rain_clusts_shorter2_no_home_dry$partid),
                      unique(Mutasa5_rain_clusts_shorter2_no_home_dry$partid),
                      unique(Nchelenge5_clusts_shorter2_no_home_dry$partid))){
    loc_count$season[ii] <- "dry"
  }
  if(temp_part %in% c(unique(Choma5_rain_clusts_shorter2_no_home_rainy$partid),
                      unique(Mutasa5_rain_clusts_shorter2_no_home_rainy$partid),
                      unique(Nchelenge5_clusts_shorter2_no_home_rainy$partid))){
    loc_count$season[ii] <- "rainy"
  }
}

loc_count$age_grp <- factor(NA, levels=c("1-16", "17-34","35-54","55+"))
loc_count$age_grp[which(loc_count$age < 17)] <- "1-16"
loc_count$age_grp[which(loc_count$age >= 17 & loc_count$age < 35)] <- "17-34"
loc_count$age_grp[which(loc_count$age >= 35 & loc_count$age < 55)] <- "35-54"
loc_count$age_grp[which(loc_count$age >= 55)] <- "55+"

loc_countb <- loc_count[which(!is.na(loc_count$age_grp)),]
loc_countb2 <- loc_countb[which(!is.na(loc_countb$male)),]
loc_countb3 <- loc_countb2[which(!is.na(loc_countb2$season)),]

loc_count2 <- loc_countb3[,c(1,5,7,3,6,2)]

descdist(loc_count2$num_locs, boot=1000, discrete = TRUE)
a1 <- fitdist(loc_count2$num_locs, "pois")
a2 <- fitdist(loc_count2$num_locs, "nbinom")

denscomp(list(a1,a2)) 
ppcomp(list(a1,a2)) ## nbinom better

bbb <- glm.nb(num_locs~loc+season+male+age_grp, data=loc_count2)
## all but season had significance

bb0 <- glm.nb(num_locs~1, data=loc_count2)
bb1 <- glm.nb(num_locs~loc, data=loc_count2)
bb2 <- glm.nb(num_locs~season, data=loc_count2)
bb3 <- glm.nb(num_locs~male, data=loc_count2)
bb4 <- glm.nb(num_locs~age_grp, data=loc_count2)
bb4b <- glm.nb(num_locs~male*age_grp, data=loc_count2)

bb5 <- glm.nb(num_locs~loc+season, data=loc_count2)
bb5b <- glm.nb(num_locs~loc*season, data=loc_count2)
bb6 <- glm.nb(num_locs~loc+male, data=loc_count2)
bb6b <- glm.nb(num_locs~loc*male, data=loc_count2)
bb7 <- glm.nb(num_locs~loc+age_grp, data=loc_count2)
bb7b <- glm.nb(num_locs~loc*age_grp, data=loc_count2)

bb8 <- glm.nb(num_locs~loc+male+age_grp, data=loc_count2)
bb8a <- glm.nb(num_locs~loc*male+age_grp, data=loc_count2)
bb8b <- glm.nb(num_locs~loc+male*age_grp, data=loc_count2)
bb8c <- glm.nb(num_locs~loc*age_grp+male, data=loc_count2)

bb9 <- glm.nb(num_locs~loc+male+season, data=loc_count2)
bb10 <- glm.nb(num_locs~loc*male*age_grp, data=loc_count2)



AIC(bb0,bb1,bb2,bb3,bb4,bb5,bb5b,bb6,bb6b,bb7,bb7b,bb8,bb8a,bb8b,bb8c,bb9,bb10)
#      df      AIC
# bb0   2 1662.562
# bb1   4 1646.127 ** #3
# bb2   3 1664.383
# bb3   3 1659.016
# bb4   5 1664.757
# bb5   5 1647.609
# bb5b  7 1649.669
# bb6   5 1642.863 ** #1
# bb6b  7 1646.078 ** #2
# bb7   7 1646.860 
# bb7b 13 1653.537

anova(bb1,bb6) ## is significantly better than jut loc to have male p=0.022
anova(bb3,bb6) ## is significantly better than just male to have loc p<0.001
### best has loc and male (but not interaction)

## is poisson better?
m3 <- glm(num_locs~loc+male, data=loc_count2, family = "poisson")
pchisq(2 * (logLik(bb6) - logLik(m3)), df = 1, lower.tail = FALSE)
## no

newdata1 <- data.frame(male = rep(c("0","1"),3), loc = rep(c("Choma","Mutasa","Nchelenge"), each=2))
newdata1$phat <- predict(bb6, newdata1, type = "response")
newdata1



dev <- c(anova(bb0, bb1, test="Chisq")[,4], anova(bb0, bb2, test="Chisq")[2,4], anova(bb0, bb3, test="Chisq")[2,4],
         anova(bb0, bb4, test="Chisq")[2,4],anova(bb0, bb4b, test="Chisq")[2,4], anova(bb0, bb5, test="Chisq")[2,4],anova(bb0, bb5b, test="Chisq")[2,4],
         anova(bb0, bb6, test="Chisq")[2,4],anova(bb0, bb6b, test="Chisq")[2,4],
         anova(bb0, bb7, test="Chisq")[2,4],anova(bb0, bb7b, test="Chisq")[2,4],
         anova(bb0, bb8, test="Chisq")[2,4],anova(bb0, bb8a, test="Chisq")[2,4],
         anova(bb0, bb8b, test="Chisq")[2,4],anova(bb0, bb8c, test="Chisq")[2,4],
         anova(bb0, bb9, test="Chisq")[2,4],anova(bb0, bb10, test="Chisq")[2,4])

b <- AICctab(bb0,bb1,bb2,bb3,bb4,bb4b,bb5,bb5b,bb6,bb6b,bb7,bb7b,bb8,bb8a,bb8b,bb8c,bb9,bb10, weights=TRUE, base=TRUE, sort=FALSE)



dev <- c(anova(bb0, bb1, test="Chisq")[,4], anova(bb0, bb2, test="Chisq")[2,4], anova(bb0, bb3, test="Chisq")[2,4],
         anova(bb0, bb4, test="Chisq")[2,4],anova(bb0, bb4b, test="Chisq")[2,4], anova(bb0, bb5, test="Chisq")[2,4],
        anova(bb0, bb6, test="Chisq")[2,4],anova(bb0, bb7, test="Chisq")[2,4],anova(bb0, bb8b, test="Chisq")[2,4])

b <- AICctab(bb0,bb1,bb2,bb3,bb4,bb4b,bb5,bb6,bb7,bb8b, weights=TRUE, base=TRUE, sort=FALSE)

cbind(data.frame(b[1:4]), dev) 



embb8 <- emmeans(bb8b, ~  loc+age_grp+male, type="response")
embb8a <- emmeans(bb8b, ~  male, type="response")
embb8b <- emmeans(bb8b, ~  loc, type="response")
embb8c <- emmeans(bb8b, ~  age_grp, type="response")

summary(bb8b)
anova(bb8b, test="Chisq")

exp(cbind(OR = coef(bb8b), confint.default(bb8b)))
# estimate the incidence of the response variable in the given category relative to the control group

emmip(bb8b, loc~male+age_grp,CIs=TRUE, type="response") +
  #scale_x_discrete(labels=c("Female","Male")) +
  #scale_y_continuous(limits=c(0,16), breaks=c(seq(0,16,2))) +
  #labs(x="Gender", y="Predicted Number of Locations Visited") +
  scale_color_manual("Site", labels=c("Choma","Mutasa","Nchelenge"),values = cols2) +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.text = element_text(size=12))




emm_dat <- emmip(bb8b, loc~male+age_grp,CIs=TRUE, type="response", plotit = FALSE) 

ggplot(data=emm_dat , aes(x=xvar, y=yvar, group=loc, color=loc)) +
  geom_point(position=position_dodge(width=0.75), size=5) +
  geom_linerange(aes(x=xvar, ymin=LCL, ymax=UCL, group=loc), 
                  position = position_dodge(width=0.75), size=7, alpha=0.4) +
  scale_x_discrete(labels=c("Female","Male","Female","Male","Female","Male","Female","Male")) +
  scale_y_continuous(limits=c(0,18.2), breaks=c(seq(0,18,2))) +
  labs(x="Gender & Age Group", y="Predicted Number of Locations Visited") +
  scale_color_manual("Site", labels=c("Choma","Mutasa","Nchelenge"),values = c(SteppedSequential5Steps[c(21,18,15)])) +
  coord_cartesian(ylim=c(0,18.2), clip="off") +
  geom_text(x=1.5,y=-2.2, label="1-16", color="black", size=5.5) +
  geom_text(x=3.5,y=-2.2, label="17-34", color="black", size=5.5) +
  geom_text(x=5.5,y=-2.2, label="35-54", color="black", size=5.5) +
  geom_text(x=7.5,y=-2.2, label="55+", color="black", size=5.5) +
  geom_vline(xintercept = 2.5, linetype=2) +
  geom_vline(xintercept = 4.5, linetype=2) +
  geom_vline(xintercept = 6.5, linetype=2) +
  theme_classic() +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=15),
        legend.title = element_text(size=14),
        axis.title.x=element_text(margin=margin(25,0,0,0)),
        legend.title.align=0.5,
        legend.text = element_text(size=12))
#

loc_count2$male_age_grp <- paste0(loc_count2$male,".",loc_count2$age_grp)
loc_count2$loc_male_age_grp <- paste0(loc_count2$loc,".",loc_count2$male,".",loc_count2$age_grp)

ggplot(data=loc_count2 , aes(x=num_locs, y=loc_male_age_grp, fill=loc, color=male)) +
  geom_density_ridges(panel_scaling = TRUE) +
  scale_fill_manual("Site", labels=c("Choma","Mutasa","Nchelenge"),values = cols2) +
  scale_color_manual("Gender", labels=c("Female","Male"),values = c("red","blue")) +
  theme_classic() +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=15),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.text = element_text(size=12))

ggplot(data=loc_count2[which(loc_count2$loc == "Mutasa" & loc_count2$male == 1),] , aes(x=num_locs, group=age_grp, color=age_grp)) +
  geom_density() +
  #scale_color_manual("Gender", labels=c("Female","Male"),values = c("red","blue")) +
  theme_classic() +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=15),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.text = element_text(size=12))

###
##
##
#############################################

medians_km_dry <- as.data.frame(cbind(c((favstats(Choma5_rain_clusts_shorter2_no_home_dry$hhdist_m_haversine~Choma5_rain_clusts_shorter2_no_home_dry$partid)[,c(1)]),
                                        (favstats(Mutasa5_rain_clusts_shorter2_no_home_dry$hhdist_m_haversine~Mutasa5_rain_clusts_shorter2_no_home_dry$partid)[,c(1)]), 
                                        (favstats(Nchelenge5_clusts_shorter2_no_home_dry$hhdist_m_haversine2~Nchelenge5_clusts_shorter2_no_home_dry$partid)[,c(1)])),
                                      c((favstats(Choma5_rain_clusts_shorter2_no_home_dry$hhdist_m_haversine~Choma5_rain_clusts_shorter2_no_home_dry$partid)[,c(4)]/1000),
                                        (favstats(Mutasa5_rain_clusts_shorter2_no_home_dry$hhdist_m_haversine~Mutasa5_rain_clusts_shorter2_no_home_dry$partid)[,c(4)]/1000), 
                                        (favstats(Nchelenge5_clusts_shorter2_no_home_dry$hhdist_m_haversine2~Nchelenge5_clusts_shorter2_no_home_dry$partid)[,c(4)]/1000)),
                                      c(rep("Choma", times=length(favstats(Choma5_rain_clusts_shorter2_no_home_dry$hhdist_m_haversine~Choma5_rain_clusts_shorter2_no_home_dry$partid)$median)),
                                        rep("Mutasa", times=length(favstats(Mutasa5_rain_clusts_shorter2_no_home_dry$hhdist_m_haversine~Mutasa5_rain_clusts_shorter2_no_home_dry$partid)$median)),
                                        rep("Nchelenge", times=length(favstats(Nchelenge5_clusts_shorter2_no_home_dry$hhdist_m_haversine~Nchelenge5_clusts_shorter2_no_home_dry$partid)$median)))))

medians_km_rainy <- as.data.frame(cbind(c((favstats(Choma5_rain_clusts_shorter2_no_home_rainy$hhdist_m_haversine~Choma5_rain_clusts_shorter2_no_home_rainy$partid)[,c(1)]),
                                          (favstats(Mutasa5_rain_clusts_shorter2_no_home_rainy$hhdist_m_haversine~Mutasa5_rain_clusts_shorter2_no_home_rainy$partid)[,c(1)]), 
                                          (favstats(Nchelenge5_clusts_shorter2_no_home_rainy$hhdist_m_haversine2~Nchelenge5_clusts_shorter2_no_home_rainy$partid)[,c(1)])),
                                        c((favstats(Choma5_rain_clusts_shorter2_no_home_rainy$hhdist_m_haversine~Choma5_rain_clusts_shorter2_no_home_rainy$partid)[,c(4)]/1000),
                                          (favstats(Mutasa5_rain_clusts_shorter2_no_home_rainy$hhdist_m_haversine~Mutasa5_rain_clusts_shorter2_no_home_rainy$partid)[,c(4)]/1000), 
                                          (favstats(Nchelenge5_clusts_shorter2_no_home_rainy$hhdist_m_haversine2~Nchelenge5_clusts_shorter2_no_home_rainy$partid)[,c(4)]/1000)),
                                        c(rep("Choma", times=length(favstats(Choma5_rain_clusts_shorter2_no_home_rainy$hhdist_m_haversine~Choma5_rain_clusts_shorter2_no_home_rainy$partid)$median)),
                                          rep("Mutasa", times=length(favstats(Mutasa5_rain_clusts_shorter2_no_home_rainy$hhdist_m_haversine~Mutasa5_rain_clusts_shorter2_no_home_rainy$partid)$median)),
                                          rep("Nchelenge", times=length(favstats(Nchelenge5_clusts_shorter2_no_home_rainy$hhdist_m_haversine~Nchelenge5_clusts_shorter2_no_home_rainy$partid)$median)))))

medians_km <- rbind(medians_km_dry, medians_km_rainy)
medians_km$season <- c(rep("dry", times=nrow(medians_km_dry)), rep("rainy", times=nrow(medians_km_rainy)))

medians_km$age <- as.numeric(NA)
medians_km$male <- factor(NA, levels=c(0,1))
for(ii in 1:length(medians_km$V1)){
  temp_part <- medians_km$V1[ii]
  a <- c(Choma5_clusts_shorter2_no_home$age_all[which(Choma5_clusts_shorter2_no_home$partid == temp_part)][1],
         Mutasa5_clusts_shorter2_no_home$age_all[which(Mutasa5_clusts_shorter2_no_home$partid == temp_part)][1],
         Nchelenge5_clusts_shorter2_no_home$age_all[which(Nchelenge5_clusts_shorter2_no_home$partid == temp_part)][1])
  a2 <- a[which(!(is.na(a)))]
  if(length(a2) > 0){
    medians_km$age[ii] <- a2
  }
  b <- c(as.character(Choma5_clusts_shorter2_no_home$male[which(Choma5_clusts_shorter2_no_home$partid == temp_part)][1]),
         as.character(Mutasa5_clusts_shorter2_no_home$male[which(Mutasa5_clusts_shorter2_no_home$partid == temp_part)][1]),
         as.character(Nchelenge5_clusts_shorter2_no_home$male[which(Nchelenge5_clusts_shorter2_no_home$partid == temp_part)][1]))
  b2 <- b[which(!(is.na(b)))]
  if(length(b2) > 0){
    medians_km$male[ii] <- b2
  }
}

medians_km$partid <- as.factor(medians_km$V1)
medians_km$dist <- as.numeric(medians_km$V2)
medians_km$loc <- as.factor(medians_km$V3)
medians_km$season <- as.factor(medians_km$season)

medians_km$age_grp <- factor(NA, levels=c("1-16", "17-34","35-54","55+"))
medians_km$age_grp[which(medians_km$age < 17)] <- "1-16"
medians_km$age_grp[which(medians_km$age >= 17 & medians_km$age < 35)] <- "17-34"
medians_km$age_grp[which(medians_km$age >= 35 & medians_km$age < 55)] <- "35-54"
medians_km$age_grp[which(medians_km$age >= 55)] <- "55+"

medians_kmb <- medians_km[which(!is.na(medians_km$age_grp)),]
medians_kmb2 <- medians_kmb[which(!is.na(medians_kmb$male)),]

medians_km2 <- medians_kmb2[,c(7,6,10,9,4,8)]
medians_km2b <- droplevels(medians_km2[which(medians_km2$dist < 85),])


descdist(medians_km2b$dist, boot=1000)
a1 <- fitdist(medians_km2b$dist, "exp")
a2 <- fitdist(medians_km2b$dist, "gamma")
a3 <- fitdist(medians_km2b$dist, "weibull")
a4 <- fitdist(medians_km2b$dist, "lnorm")
a5 <- fitdist(medians_km2b$dist, "norm")

denscomp(list(a1,a2,a3,a4,a5)) 
ppcomp(list(a1,a2,a3,a4,a5)) ## lnorm better

summary(a4)
#          estimate Std. Error
# meanlog 0.4689863 0.08332507
# sdlog   1.3487341 0.05891958


##
a <- gamlss(dist~1, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO)

a1 <- GAIC(gamlss(dist~1, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~season, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~male, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO))
# + loc ## AIC 16 less, df 1 higher

a2 <- GAIC(gamlss(dist~loc, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+season, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+male, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO))
# + loc + season ## AIC 6 less, df 1 higher

a3 <- GAIC(gamlss(dist~loc+season, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+season+male, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+season+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO))
# + loc + season + age_grp ## AIC down by 6.7 df up by 3

a4 <- GAIC(gamlss(dist~loc+season+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+season+age_grp+male, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO))
# + loc + season + age_grp 

## SIGMA STEP
a5 <- GAIC(gamlss(dist~loc+season+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+season+age_grp, sigma.formula = ~ loc, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+season+age_grp, sigma.formula = ~ season, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+season+age_grp, sigma.formula = ~ male, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+season+age_grp, sigma.formula = ~ age_grp, data=medians_km2b, family=LOGNO))
# sigma = ~age_grp ## AIC down by 0.4 df up by 3 --> not worth it

## MU STEP DOWN
a6 <- GAIC(gamlss(dist~loc+season+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+season, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~loc+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO),
           gamlss(dist~season+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO))
# + loc + season + age_grp ## AIC down by 6.7 df up by 3 --> still on the fence


aa10 <- gamlss(dist~loc+season+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO)
aa11 <- gamlss(dist~loc+season, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO)
aa12 <- gamlss(dist~loc+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO)
aa13 <- gamlss(dist~season+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO)
aa10a <- gamlss(dist~loc*season+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO)
aa10b <- gamlss(dist~loc+season*age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO)
aa10c <- gamlss(dist~loc*age_grp+season, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO)
GAIC(aa10,aa11,aa12,aa13,aa10a,aa10b,aa10c)
## a10a has AIC 2 lower but df 2 higher --> not better


# aa10 <- gamlss(dist~loc+season+age_grp, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO)
newdata2 <- data.frame(season = rep(c(rep(c("dry","rainy"),each=4)),3), 
                       loc = rep(c("Choma","Mutasa","Nchelenge"), each=8), 
                       age_grp = rep(c("1-16","17-34","35-54","55+"), 6))
q <- predictAll(aa10, newdata = newdata2, type="response", se.fit = TRUE) # fitted 
q$fitted_exp_val <- exp(q$mu + (q$sigma^2)/2)
q$fitted_sd_val <- sqrt((exp(q$sigma^2) - 1)*(exp(2*q$mu + q$sigma^2)))

q$fitted_exp_val_low_sd <- q$fitted_exp_val - q$fitted_sd_val
q$fitted_exp_val_high_sd <- q$fitted_exp_val + q$fitted_sd_val

q$mu_low <- q$mu - q$mu.se
q$mu_high <- q$mu + q$mu.se

q$fitted_exp_val_low <- exp(q$mu_low + (q$sigma^2)/2)
q$fitted_exp_val_high <- exp(q$mu_high + (q$sigma^2)/2)



q2 <- cbind(newdata2,q)
q2$season <- as.factor(q2$season)
q2$loc <- as.factor(q2$loc)
q2$age_grp <- as.factor(q2$age_grp)
q2$loc_season <- factor("Choma_dry", levels=c("Choma_dry", "Choma_rainy", "Mutasa_dry","Mutasa_rainy",
                                        "Nchelenge_dry", "Nchelenge_rainy"))
q2$loc_season[which(q2$loc == "Choma" & q2$season == "rainy")] <- "Choma_rainy"
q2$loc_season[which(q2$loc == "Mutasa" & q2$season == "dry")] <- "Mutasa_dry"
q2$loc_season[which(q2$loc == "Mutasa" & q2$season == "rainy")] <- "Mutasa_rainy"
q2$loc_season[which(q2$loc == "Nchelenge" & q2$season == "dry")] <- "Nchelenge_dry"
q2$loc_season[which(q2$loc == "Nchelenge" & q2$season == "rainy")] <- "Nchelenge_rainy"



ggplot(data=q2 , aes(x=loc_season, y=fitted_exp_val, group=age_grp, color=age_grp)) +
  
  geom_linerange(aes(x=loc_season, ymin=fitted_exp_val_low, ymax=fitted_exp_val_high, group=age_grp, color=age_grp), 
                 position = position_dodge(width=0.75), size=7, alpha=0.6) +
  
  geom_point(size=5,position = position_dodge(width=0.75)) +
  
  
  scale_x_discrete(labels=c("dry", "rainy","dry", "rainy","dry", "rainy")) +
  scale_y_continuous(limits=c(0,20), breaks=c(seq(0,20,2))) +
  labs(y="Predicted Median Distance of Locations Visited", x="Site and Season") +
  scale_color_manual("Age", labels=c("1-16","17-34","35-54","55+"),values = agecols) +
  coord_cartesian(ylim=c(0,17.5), clip="off") +
  geom_text(x=1.5,y=-2, label="Choma", color="black", size=5) +
  geom_text(x=3.5,y=-2, label="Mutasa", color="black", size=5) +
  geom_text(x=5.5,y=-2, label="Nchelenge", color="black", size=5) +
  geom_vline(xintercept = 2.5, linetype=2) +
  geom_vline(xintercept = 4.5, linetype=2) +
  theme_classic() +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=17),
        axis.title.x=element_text(margin=margin(20,0,0,0)),
        legend.title = element_text(size=16),
        legend.title.align=0.5,
        legend.text = element_text(size=14))




## next best
# aa11 <- gamlss(dist~loc+season, sigma.formula = ~ 1, data=medians_km2b, family=LOGNO)
newdata2 <- data.frame(season = rep(c("dry","rainy"),3), 
                       loc = rep(c("Choma","Mutasa","Nchelenge"), each=2))
q <- predictAll(aa11, newdata = newdata2, type="response", se.fit = TRUE) # fitted 
q$fitted_exp_val <- exp(q$mu + (q$sigma^2)/2)
q$fitted_sd_val <- sqrt((exp(q$sigma^2) - 1)*(exp(2*q$mu + q$sigma^2)))

q$fitted_exp_val_low_sd <- q$fitted_exp_val - q$fitted_sd_val
q$fitted_exp_val_high_sd <- q$fitted_exp_val + q$fitted_sd_val

q$mu_low <- q$mu - q$mu.se
q$mu_high <- q$mu + q$mu.se

q$fitted_exp_val_low <- exp(q$mu_low + (q$sigma^2)/2)
q$fitted_exp_val_high <- exp(q$mu_high + (q$sigma^2)/2)




q2 <- cbind(newdata2,q)
q2$season <- as.factor(q2$season)
q2$loc <- as.factor(q2$loc)

ggplot(data=q2 , aes(x=season, y=fitted_exp_val, group=loc, color=loc)) +
  geom_point(size=5,position = position_dodge(width=0.75)) +
  
  geom_linerange(aes(x=season, ymin=fitted_exp_val_low, ymax=fitted_exp_val_high, group=loc, color=loc), 
                 position = position_dodge(width=0.75), size=7, alpha=0.4) +
  
  scale_x_discrete(labels=c("Dry", "Rainy")) +
  scale_y_continuous(limits=c(0,16), breaks=c(seq(0,20,2))) +
  labs(y="Predicted Distance of Locations Visited", x="Season") +
  scale_color_manual("Site", labels=c("Choma","Mutasa","Nchelenge"),values = cols2) +
  theme_classic() +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=15),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.text = element_text(size=12))

### Mutasa and Nchelenge just slightly different ####
########
########
#### COMPOSITE VARIABLE ####
## score per individual based on percentile/quartile fall into for # locations, # trips, distance (max or average)
# add up percentiles for different individuals to make summary metric where larger values are those people who are highly mobile
# basic regression to predict value based on age, gender, site to see if general patterns of who is highly mobile or not

# loc_count2 has gender, age, loc, season and # locs --> add # trips, distance
loc_counts_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust~Choma5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,9)]
loc_counts_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,9)]

only_home_Choma <- levels(as.factor(Choma5_clusts_shorter2$partid))[which(!(levels(as.factor(Choma5_clusts_shorter2$partid)) %in% levels(as.factor(Choma5_clusts_shorter2_no_home$partid))))]
only_home_Mutasa <- levels(as.factor(Mutasa5_clusts_shorter2$partid))[which(!(levels(as.factor(Mutasa5_clusts_shorter2$partid)) %in% levels(as.factor(Mutasa5_clusts_shorter2_no_home$partid))))]
only_home_Nchelenge <- levels(as.factor(Nchelenge5_clusts_shorter2$partid))[which(!(levels(as.factor(Nchelenge5_clusts_shorter2$partid)) %in% levels(as.factor(Nchelenge5_clusts_shorter2_no_home$partid))))]


loc_count <- as.data.frame(cbind(c(loc_counts_Choma$`Choma5_clusts_shorter2_no_home$partid`, 
                                   loc_counts_Mutasa$`Mutasa5_clusts_shorter2_no_home$partid`, 
                                   loc_counts_Nchelenge$`Nchelenge5_clusts_shorter2_no_home$partid`,
                                   only_home_Choma, only_home_Mutasa, only_home_Nchelenge),
                                 c(loc_counts_Choma$n, 
                                   loc_counts_Mutasa$n, 
                                   loc_counts_Nchelenge$n,
                                   rep(0,length(only_home_Choma)+length(only_home_Mutasa)+length(only_home_Nchelenge))),
                                 c(rep("Choma",length(loc_counts_Choma$n)),
                                   rep("Mutasa",length(loc_counts_Mutasa$n)),
                                   rep("Nchelenge",length(loc_counts_Nchelenge$n)),
                                   rep("Choma",length(only_home_Choma)),
                                   rep("Mutasa",length(only_home_Mutasa)),
                                   rep("Nchelenge",length(only_home_Nchelenge)))))
colnames(loc_count) <- c("partid","num_locs","loc")
loc_count$num_locs <- as.numeric(as.character(loc_count$num_locs))
loc_count$partid <- factor(loc_count$partid)
loc_count$loc <- factor(loc_count$loc)

loc_count$age <- as.numeric(NA)
loc_count$male <- factor(NA, levels=c(0,1))
loc_count$season <- factor(NA, levels=c("dry","rainy"))
for(ii in 1:length(loc_count$partid)){
  temp_part <- loc_count$partid[ii]
  a <- c(Choma5_clusts_shorter2$age_all[which(Choma5_clusts_shorter2$partid == temp_part)][1],
         Mutasa5_clusts_shorter2$age_all[which(Mutasa5_clusts_shorter2$partid == temp_part)][1],
         Nchelenge5_clusts_shorter2$age_all[which(Nchelenge5_clusts_shorter2$partid == temp_part)][1])
  a2 <- a[which(!(is.na(a)))]
  if(length(a2) > 0){
    loc_count$age[ii] <- a2
  }
  b <- c(as.character(Choma5_clusts_shorter2$male[which(Choma5_clusts_shorter2$partid == temp_part)][1]),
         as.character(Mutasa5_clusts_shorter2$male[which(Mutasa5_clusts_shorter2$partid == temp_part)][1]),
         as.character(Nchelenge5_clusts_shorter2$male[which(Nchelenge5_clusts_shorter2$partid == temp_part)][1]))
  b2 <- b[which(!(is.na(b)))]
  if(length(b2) > 0){
    loc_count$male[ii] <- b2
  }
  if(temp_part %in% c(unique(Choma5_rain_clusts_shorter2_dry$partid),
                      unique(Mutasa5_rain_clusts_shorter2_dry$partid),
                      unique(Nchelenge5_clusts_shorter2_dry$partid))){
    loc_count$season[ii] <- "dry"
  }
  if(temp_part %in% c(unique(Choma5_rain_clusts_shorter2_rainy$partid),
                      unique(Mutasa5_rain_clusts_shorter2_rainy$partid),
                      unique(Nchelenge5_clusts_shorter2_rainy$partid))){
    loc_count$season[ii] <- "rainy"
  }
}

loc_count$age_grp <- factor(NA, levels=c("1-16", "17-34","35-54","55+"))
loc_count$age_grp[which(loc_count$age < 17)] <- "1-16"
loc_count$age_grp[which(loc_count$age >= 17 & loc_count$age < 35)] <- "17-34"
loc_count$age_grp[which(loc_count$age >= 35 & loc_count$age < 55)] <- "35-54"
loc_count$age_grp[which(loc_count$age >= 55)] <- "55+"

loc_countb <- loc_count[which(!is.na(loc_count$age_grp)),]
loc_countb2 <- loc_countb[which(!is.na(loc_countb$male)),]
loc_countb3 <- loc_countb2[which(!is.na(loc_countb2$season)),]

loc_count2 <- loc_countb3[,c(1,5,7,3,6,2)]

loc_all <- loc_count2

trips_part_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home3$trip_count_new ~ Choma5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]
trips_part_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home3$trip_count_new ~ Mutasa5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]
trips_part_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home3$trip_count_new ~ Nchelenge5_clust_trips_shorter2_no_home3$partid)[,c(1,4)]
names(trips_part_median_Choma) <- c("partid", "median_trips")
names(trips_part_median_Mutasa) <- c("partid", "median_trips")
names(trips_part_median_Nchelenge) <- c("partid", "median_trips")
trips_part_medians <- rbind(trips_part_median_Choma,trips_part_median_Mutasa,trips_part_median_Nchelenge)
loc_all2 <- merge(loc_all, trips_part_medians, by="partid", all.x=TRUE)


medians_km_Choma <- favstats(Choma5_clusts_shorter2_no_home$hhdist_km_haversine~Choma5_clusts_shorter2_no_home$partid)[,c(1,4)]
medians_km_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,4)]
medians_km_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,4)]
names(medians_km_Choma) <- c("partid", "median_dist")
names(medians_km_Mutasa) <- c("partid", "median_dist")
names(medians_km_Nchelenge) <- c("partid", "median_dist")
dist_part_medians <- rbind(medians_km_Choma,medians_km_Mutasa,medians_km_Nchelenge)
loc_all3 <- merge(loc_all2, dist_part_medians, by="partid",all.x=TRUE)


max_km_Choma <- favstats(Choma5_clusts_shorter2_no_home$hhdist_km_haversine~Choma5_clusts_shorter2_no_home$partid)[,c(1,6)]
max_km_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$hhdist_km_haversine~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,6)]
max_km_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$hhdist_km_haversine~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,6)]
names(max_km_Choma) <- c("partid", "max_dist")
names(max_km_Mutasa) <- c("partid", "max_dist")
names(max_km_Nchelenge) <- c("partid", "max_dist")
dist_part_max <- rbind(max_km_Choma,max_km_Mutasa,max_km_Nchelenge)
loc_all4 <- merge(loc_all3, dist_part_max, by="partid",all.x=TRUE)

clusts_time_median_Choma <- favstats(Choma5_clusts_shorter2_no_home$clust_trips_time_new~Choma5_clusts_shorter2_no_home$partid)[,c(1,4)]
clusts_time_median_Mutasa <- favstats(Mutasa5_clusts_shorter2_no_home$clust_trips_time_new~Mutasa5_clusts_shorter2_no_home$partid)[,c(1,4)]
clusts_time_median_Nchelenge <- favstats(Nchelenge5_clusts_shorter2_no_home$clust_trips_time_new~Nchelenge5_clusts_shorter2_no_home$partid)[,c(1,4)]
names(clusts_time_median_Choma) <- c("partid", "med_time")
names(clusts_time_median_Mutasa) <- c("partid", "med_time")
names(clusts_time_median_Nchelenge) <- c("partid", "med_time")
clust_time_part_med <- rbind(clusts_time_median_Choma,clusts_time_median_Mutasa,clusts_time_median_Nchelenge)
clust_time_part_med$med_time_hr <- clust_time_part_med$med_time*24
clust_time_part_med <- clust_time_part_med[,c(1,3)]
loc_all5 <- merge(loc_all4, clust_time_part_med, by="partid",all.x=TRUE)


trip_time_median_Choma <- favstats(Choma5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Choma5_clust_trips_shorter2_no_home$partid)[,c(1,4)]
trip_time_median_Mutasa <- favstats(Mutasa5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Mutasa5_clust_trips_shorter2_no_home$partid)[,c(1,4)]
trip_time_median_Nchelenge <- favstats(Nchelenge5_clust_trips_shorter2_no_home$trip_time_new_hour ~ Nchelenge5_clust_trips_shorter2_no_home$partid)[,c(1,4)]
names(trip_time_median_Choma) <- c("partid", "med_trip_time")
names(trip_time_median_Mutasa) <- c("partid", "med_trip_time")
names(trip_time_median_Nchelenge) <- c("partid", "med_trip_time")
trip_time_part_med <- rbind(trip_time_median_Choma,trip_time_median_Mutasa,trip_time_median_Nchelenge)
trip_time_part_med$med_trip_time_hr <- trip_time_part_med$med_trip_time*24
trip_time_part_med <- trip_time_part_med[,c(1,3)]
loc_all6 <- merge(loc_all5, trip_time_part_med, by="partid",all.x=TRUE)


perc_home_Choma <- Choma5_clusts_shorter2_home[,c(1,40)]
perc_home_Mutasa <- Mutasa5_clusts_shorter2_home[,c(1,41)]
perc_home_Nchelenge <- Nchelenge5_clusts_shorter2_home[,c(1,56)]
perc_home_all <- rbind(perc_home_Choma,perc_home_Mutasa,perc_home_Nchelenge)
loc_all7 <- merge(loc_all6, perc_home_all, by="partid",all.x=TRUE)
loc_all7$perc_out <- 100-loc_all7$perc_time_home

loc_all7$quart_locs <- ntiles(loc_all7$num_locs, n=4)
loc_all7$quart_trips <- ntiles(loc_all7$median_trips, n=4)
loc_all7$quart_med_dist <- ntiles(loc_all7$median_dist, n=4)
loc_all7$quart_max_dist <- ntiles(loc_all7$max_dist, n=4)
loc_all7$quart_time <- ntiles(loc_all7$med_time_hr, n=4)
loc_all7$quart_trip_time <- ntiles(loc_all7$med_trip_time_hr, n=4)
loc_all7$quart_perc_out <- ntiles(loc_all7$perc_out, n=4)

loc_all7$quart_added_med_dist_clust_time <- as.numeric(loc_all7$quart_locs)  + as.numeric(loc_all7$quart_time) + as.numeric(loc_all7$quart_med_dist) ## of 12
loc_all7$quart_added_max_dist_clust_time <- as.numeric(loc_all7$quart_locs) + as.numeric(loc_all7$quart_time) + as.numeric(loc_all7$quart_max_dist) ## of 12
loc_all7$quart_added_both_dist_clust_time <- as.numeric(loc_all7$quart_locs) + as.numeric(loc_all7$quart_time) + as.numeric(loc_all7$quart_med_dist) + as.numeric(loc_all7$quart_max_dist) ## of 16

loc_all7$quart_added_med_dist_trips <- as.numeric(loc_all7$quart_locs) + as.numeric(loc_all7$quart_trips)  + as.numeric(loc_all7$quart_trip_time) + as.numeric(loc_all7$quart_med_dist) ## of 16
loc_all7$quart_added_max_dist_trips <- as.numeric(loc_all7$quart_locs) + as.numeric(loc_all7$quart_trips) + as.numeric(loc_all7$quart_trip_time) + as.numeric(loc_all7$quart_max_dist) ## of 16
loc_all7$quart_added_both_dist_trips <- as.numeric(loc_all7$quart_locs) + as.numeric(loc_all7$quart_trips) + as.numeric(loc_all7$quart_trip_time) + as.numeric(loc_all7$quart_med_dist) + as.numeric(loc_all7$quart_max_dist) ## of 20

loc_all7$quart_added_med_dist_all <- as.numeric(loc_all7$quart_locs)  + as.numeric(loc_all7$quart_time) + as.numeric(loc_all7$quart_trips)  + as.numeric(loc_all7$quart_trip_time) + as.numeric(loc_all7$quart_med_dist) ## of 20
loc_all7$quart_added_max_dist_all <- as.numeric(loc_all7$quart_locs) + as.numeric(loc_all7$quart_time) + as.numeric(loc_all7$quart_trips)  + as.numeric(loc_all7$quart_trip_time) + as.numeric(loc_all7$quart_max_dist) ## of 20
loc_all7$quart_added_both_dist_all <- as.numeric(loc_all7$quart_locs) + as.numeric(loc_all7$quart_time) + as.numeric(loc_all7$quart_trips)  + as.numeric(loc_all7$quart_trip_time) + as.numeric(loc_all7$quart_med_dist) + as.numeric(loc_all7$quart_max_dist) ## of 24

loc_all7$quart_added_med_dist_most <- as.numeric(loc_all7$quart_locs)  + as.numeric(loc_all7$quart_time) + as.numeric(loc_all7$quart_trips)  + as.numeric(loc_all7$quart_med_dist) ## of 16
loc_all7$quart_added_med_dist_most_perc <- as.numeric(loc_all7$quart_locs)  + as.numeric(loc_all7$quart_time) + as.numeric(loc_all7$quart_trips)  + as.numeric(loc_all7$quart_med_dist) + as.numeric(loc_all7$quart_perc_out) ## of 20
loc_all7$quart_added_med_dist_most2 <- loc_all7$quart_added_med_dist_most-4
loc_all7$quart_added_med_dist_most_perc2 <- loc_all7$quart_added_med_dist_most_perc-5
## % time outside house very positively correlated with # locs

ggplot(data=loc_all7, aes(x=quart_added_med_dist_clust_time)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,12,1))) 
ggplot(data=loc_all7, aes(x=quart_added_med_dist_trips)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,16,1))) 
ggplot(data=loc_all7, aes(x=quart_added_med_dist_all)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,20,1))) 

ggplot(data=loc_all7, aes(x=quart_added_max_dist_clust_time)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,12,1))) 
ggplot(data=loc_all7, aes(x=quart_added_max_dist_trips)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,16,1))) 
ggplot(data=loc_all7, aes(x=quart_added_max_dist_all)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,20,1))) 

ggplot(data=loc_all7, aes(x=quart_added_both_dist_clust_time)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,16,1))) 
ggplot(data=loc_all7, aes(x=quart_added_both_dist_trips)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,20,1))) 
ggplot(data=loc_all7, aes(x=quart_added_both_dist_all)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,24,1))) 

ggplot(data=loc_all7, aes(x=quart_added_med_dist_most)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,16,1))) 
ggplot(data=loc_all7, aes(x=quart_added_med_dist_most_perc2)) +
  geom_histogram(binwidth=1) +
  scale_x_continuous(breaks=c(seq(0,20,1))) 




ggplot(data=loc_all7, aes(x=quart_added_med_dist_most,  y=..ndensity.., group=male, color=male)) +
  #geom_histogram(binwidth=1, position = "dodge") +
  geom_density() +
  theme_classic() +
  scale_x_continuous(breaks=c(seq(0,20,1))) 

ggplot(data=loc_all7, aes(x=quart_added_med_dist_most,  y=..ndensity.., group=age_grp, color=age_grp)) +
  #geom_histogram(binwidth=1, position = "dodge") +
  geom_density() +
  theme_classic() +
  scale_x_continuous(breaks=c(seq(0,20,1))) 

ggplot(data=loc_all7, aes(x=quart_added_med_dist_most,  y=..ndensity.., group=loc, color=loc)) +
  #geom_histogram(binwidth=1, position = "dodge") +
  geom_density() +
  theme_classic() +
  scale_x_continuous(breaks=c(seq(0,20,1))) 

ggplot(data=loc_all7, aes(x=quart_added_med_dist_most, y=..ndensity.., group=season, color=season)) +
  #geom_histogram(binwidth=1, position = "dodge") +
  geom_density() +
  theme_classic() +
  scale_x_continuous(breaks=c(seq(0,20,1))) 


loc_all71 <- loc_all7[which(!(is.na(loc_all7$quart_added_med_dist_most_perc))),]

shapiro.test(loc_all71$quart_added_med_dist_all) ## not normal
descdist(loc_all71$quart_added_med_dist_all, boot = 1000)

a1 <- fitdist(loc_all71$quart_added_med_dist_all, "norm")
a2 <- fitdist(loc_all71$quart_added_med_dist_all, "nbinom")
a3 <- fitdist(loc_all71$quart_added_med_dist_all, "pois")
a4 <- fitdist(loc_all71$quart_added_med_dist_all, "lnorm")
a5 <- fitdist(loc_all71$quart_added_med_dist_all, "gamma")
a6 <- fitdist(loc_all71$quart_added_med_dist_all, "unif")

denscomp(list(a1,a2,a3,a4,a5,a6)) 
ppcomp(list(a1,a2,a3,a4,a5,a6)) #poisson

m_all <- glm(quart_added_med_dist_all~male+age_grp+loc+season, data=loc_all71, family=poisson)
## maybe male and maybe loc
m0 <- glm(quart_added_med_dist_all~1, data=loc_all71, family=poisson)
m1 <- glm(quart_added_med_dist_all~male, data=loc_all71, family=poisson)
m2 <- glm(quart_added_med_dist_all~age_grp, data=loc_all71, family=poisson)
m3 <- glm(quart_added_med_dist_all~loc, data=loc_all71, family=poisson)
m4 <- glm(quart_added_med_dist_all~season, data=loc_all71, family=poisson)
AIC(m0,m1,m2,m3,m4)

m5 <- glm(quart_added_med_dist_all~loc+male, data=loc_all71, family=poisson)
m6 <- glm(quart_added_med_dist_all~loc+age_grp, data=loc_all71, family=poisson)
m7 <- glm(quart_added_med_dist_all~loc+season, data=loc_all71, family=poisson)
AIC(m0,m1,m2,m3,m4,m5,m6,m7)

m8 <- glm(quart_added_med_dist_all~loc*male, data=loc_all71, family=poisson)
m9 <- glm(quart_added_med_dist_all~loc*age_grp, data=loc_all71, family=poisson)
m10 <- glm(quart_added_med_dist_all~loc*season, data=loc_all71, family=poisson)
AIC(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
## loc
b <- AICctab(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10, weights=TRUE, base=TRUE, sort=FALSE)


newdata1 <- data.frame(loc = c("Choma","Mutasa","Nchelenge"))
newdata1$phat <- predict(m3, newdata1, type = "response")
newdata1

emmip(m3, ~loc,CIs=TRUE, type="response") +
  scale_x_discrete(labels=c("Choma","Mutasa","Nchelenge")) +
  scale_y_continuous(limits=c(0,20), breaks=c(seq(0,20,2))) +
  labs(x="Site", y="Predicted Composite Score (of 20)") +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.text = element_text(size=12))








m_all <- glm(quart_added_med_dist_clust_time~male+age_grp+loc+season, data=loc_all71, family=poisson)
## maybe male and maybe loc
m0 <- glm(quart_added_med_dist_clust_time~1, data=loc_all71, family=poisson)
m1 <- glm(quart_added_med_dist_clust_time~male, data=loc_all71, family=poisson)
m2 <- glm(quart_added_med_dist_clust_time~age_grp, data=loc_all71, family=poisson)
m3 <- glm(quart_added_med_dist_clust_time~loc, data=loc_all71, family=poisson)
m4 <- glm(quart_added_med_dist_clust_time~season, data=loc_all71, family=poisson)
AIC(m0,m1,m2,m3,m4)

m5 <- glm(quart_added_med_dist_clust_time~loc+male, data=loc_all71, family=poisson)
m6 <- glm(quart_added_med_dist_clust_time~loc+age_grp, data=loc_all71, family=poisson)
m7 <- glm(quart_added_med_dist_clust_time~loc+season, data=loc_all71, family=poisson)
AIC(m0,m1,m2,m3,m4,m5,m6,m7)

m8 <- glm(quart_added_med_dist_clust_time~loc*male, data=loc_all71, family=poisson)
m9 <- glm(quart_added_med_dist_clust_time~loc*age_grp, data=loc_all71, family=poisson)
m10 <- glm(quart_added_med_dist_clust_time~loc*season, data=loc_all71, family=poisson)
AIC(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
## loc
b <- AICctab(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10, weights=TRUE, base=TRUE, sort=FALSE)


newdata1 <- data.frame(loc = c("Choma","Mutasa","Nchelenge"))
newdata1$phat <- predict(m3, newdata1, type = "response")
newdata1

emmip(m3, ~loc,CIs=TRUE, type="response") +
  scale_x_discrete(labels=c("Choma","Mutasa","Nchelenge")) +
  scale_y_continuous(limits=c(0,12), breaks=c(seq(0,12,2))) +
  labs(x="Site", y="Predicted Composite Score (of 12)") +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.text = element_text(size=12))











m_all <- glm(quart_added_med_dist_trips~male+age_grp+loc+season, data=loc_all71, family=poisson)
## maybe male and maybe loc
m0 <- glm(quart_added_med_dist_trips~1, data=loc_all71, family=poisson)
m1 <- glm(quart_added_med_dist_trips~male, data=loc_all71, family=poisson)
m2 <- glm(quart_added_med_dist_trips~age_grp, data=loc_all71, family=poisson)
m3 <- glm(quart_added_med_dist_trips~loc, data=loc_all71, family=poisson)
m4 <- glm(quart_added_med_dist_trips~season, data=loc_all71, family=poisson)
AIC(m0,m1,m2,m3,m4)

m5 <- glm(quart_added_med_dist_trips~loc+male, data=loc_all71, family=poisson)
m6 <- glm(quart_added_med_dist_trips~loc+age_grp, data=loc_all71, family=poisson)
m7 <- glm(quart_added_med_dist_trips~loc+season, data=loc_all71, family=poisson)
AIC(m0,m1,m2,m3,m4,m5,m6,m7)

m8 <- glm(quart_added_med_dist_trips~loc*male, data=loc_all71, family=poisson)
m9 <- glm(quart_added_med_dist_trips~loc*age_grp, data=loc_all71, family=poisson)
m10 <- glm(quart_added_med_dist_trips~loc*season, data=loc_all71, family=poisson)
AIC(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
## loc
b <- AICctab(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10, weights=TRUE, base=TRUE, sort=FALSE)


newdata1 <- data.frame(loc = c("Choma","Mutasa","Nchelenge"))
newdata1$phat <- predict(m3, newdata1, type = "response")
newdata1

emmip(m3, ~loc,CIs=TRUE, type="response") +
  scale_x_discrete(labels=c("Choma","Mutasa","Nchelenge")) +
  scale_y_continuous(limits=c(0,16), breaks=c(seq(0,16,2))) +
  labs(x="Site", y="Predicted Composite Score (of 16)") +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.text = element_text(size=12))







###### option 1 ####

# m_all <- glm(quart_added_med_dist_most2~male+age_grp+loc+season, data=loc_all71, family=poisson)
# ## maybe male and maybe loc
# m0 <- glm(quart_added_med_dist_most2~1, data=loc_all71, family=poisson)
# m1 <- glm(quart_added_med_dist_most2~male, data=loc_all71, family=poisson)
# m2 <- glm(quart_added_med_dist_most2~age_grp, data=loc_all71, family=poisson)
# m3 <- glm(quart_added_med_dist_most2~loc, data=loc_all71, family=poisson)
# m4 <- glm(quart_added_med_dist_most2~season, data=loc_all71, family=poisson)
# AIC(m0,m1,m2,m3,m4)
# 
# m5 <- glm(quart_added_med_dist_most2~loc+male, data=loc_all71, family=poisson)
# m6 <- glm(quart_added_med_dist_most2~loc+age_grp, data=loc_all71, family=poisson)
# m7 <- glm(quart_added_med_dist_most2~loc+season, data=loc_all71, family=poisson)
# AIC(m0,m1,m2,m3,m4,m5,m6,m7)
# 
# m8 <- glm(quart_added_med_dist_most2~loc*male, data=loc_all71, family=poisson)
# m9 <- glm(quart_added_med_dist_most2~loc*age_grp, data=loc_all71, family=poisson)
# m10 <- glm(quart_added_med_dist_most2~loc*season, data=loc_all71, family=poisson)
# AIC(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
# ## loc
# b <- AICctab(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10, weights=TRUE, base=TRUE, sort=FALSE)
# ## loc+male then loc*male then loc
# 
# newdata1 <- data.frame(loc = rep(c("Choma","Mutasa","Nchelenge"),2), male=rep(c("0","1"),each=3))
# newdata1$phat <- predict(m5, newdata1, type = "response")
# newdata1
# 
# emmip(m5, loc~male,CIs=TRUE, type="response") +
#   scale_x_discrete(labels=c("Female","Male")) +
#   scale_y_continuous(limits=c(0,12), breaks=c(seq(0,12,2))) +
#   labs(x="Gender", y="Predicted Composite Score (of 12)") +
#   scale_color_manual("Site", labels=c("Choma","Mutasa","Nchelenge"),values = cols2) +
#   theme_bw() +
#   theme(axis.text = element_text(size=14),
#         axis.title = element_text(size=15),
#         legend.title = element_text(size=14),
#         legend.title.align=0.5,
#         legend.text = element_text(size=12))
# 
# emmip(m8, loc~male,CIs=TRUE, type="response") +
#   scale_x_discrete(labels=c("Female","Male")) +
#   scale_y_continuous(limits=c(0,12), breaks=c(seq(0,12,2))) +
#   labs(x="Gender", y="Predicted Composite Score (of 12)") +
#   scale_color_manual("Site", labels=c("Choma","Mutasa","Nchelenge"),values = cols2) +
#   theme_bw() +
#   theme(axis.text = element_text(size=14),
#         axis.title = element_text(size=15),
#         legend.title = element_text(size=14),
#         legend.title.align=0.5,
#         legend.text = element_text(size=12))
# 
# emmip(m3, ~loc,CIs=TRUE, type="response") +
#   scale_x_discrete(labels=c("Choma","Mutasa","Nchelenge")) +
#   scale_y_continuous(limits=c(0,12), breaks=c(seq(0,12,2))) +
#   labs(x="Site", y="Predicted Composite Score (of 12)") +
#   theme_bw() +
#   theme(axis.text = element_text(size=14),
#         axis.title = element_text(size=15),
#         legend.title = element_text(size=14),
#         legend.title.align=0.5,
#         legend.text = element_text(size=12))
# 
# 

##### option 2 ####

m_all <- glm(quart_added_med_dist_most_perc2~male+age_grp+loc+season, data=loc_all71, family=poisson)
## maybe male and maybe loc
m0 <- glm(quart_added_med_dist_most_perc2~1, data=loc_all71, family=poisson)
m1 <- glm(quart_added_med_dist_most_perc2~male, data=loc_all71, family=poisson)
m2 <- glm(quart_added_med_dist_most_perc2~age_grp, data=loc_all71, family=poisson)
m3 <- glm(quart_added_med_dist_most_perc2~loc, data=loc_all71, family=poisson)
m4 <- glm(quart_added_med_dist_most_perc2~season, data=loc_all71, family=poisson)
AIC(m0,m1,m2,m3,m4)

m5 <- glm(quart_added_med_dist_most_perc2~loc+male, data=loc_all71, family=poisson)
m6 <- glm(quart_added_med_dist_most_perc2~loc+age_grp, data=loc_all71, family=poisson)
m7 <- glm(quart_added_med_dist_most_perc2~loc+season, data=loc_all71, family=poisson)
AIC(m0,m1,m2,m3,m4,m5,m6,m7)

m8 <- glm(quart_added_med_dist_most_perc2~loc*male, data=loc_all71, family=poisson)
m9 <- glm(quart_added_med_dist_most_perc2~loc*age_grp, data=loc_all71, family=poisson)
m10 <- glm(quart_added_med_dist_most_perc2~loc*season, data=loc_all71, family=poisson)

AIC(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
## loc
b <- AICctab(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10, weights=TRUE, base=TRUE, sort=FALSE)
## loc*male then loc+male then loc*season

newdata1 <- data.frame(loc = rep(c("Choma","Mutasa","Nchelenge"),2), male=rep(c("0","1"),each=3))
newdata1$phat <- predict(m8, newdata1, type = "response")
newdata1
emmeans(m8, ~  loc, type="response")
emmip(m8, loc~male,CIs=TRUE, type="response") +
  scale_x_discrete(labels=c("Female","Male")) +
  scale_y_continuous(limits=c(0,15), breaks=c(seq(0,15,2))) +
  labs(x="Gender", y="Predicted Composite Score (of 15)") +
  scale_color_manual("Site", labels=c("Choma","Mutasa","Nchelenge"),values = cols2) +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.text = element_text(size=12))

emmip(m5, loc~male,CIs=TRUE, type="response") +
  scale_x_discrete(labels=c("Female","Male")) +
  scale_y_continuous(limits=c(0,15), breaks=c(seq(0,15,2))) +
  labs(x="Gender", y="Predicted Composite Score (of 15)") +
  scale_color_manual("Site", labels=c("Choma","Mutasa","Nchelenge"),values = cols2) +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.text = element_text(size=12))
##


emm_dat <- emmip(m8, loc~male,CIs=TRUE, type="response", plotit = FALSE) 

ggplot(data=emm_dat , aes(x=yvar, y=xvar, group=loc, color=loc)) +
  geom_point(position=position_dodge(width=0.75), size=5) +
  geom_linerange(aes(x=yvar, xmin=LCL, xmax=UCL, group=loc), 
                 position = position_dodge(width=0.75), size=7, alpha=0.4) +
  scale_y_discrete(labels=c("Female","Male","Female","Male","Female","Male","Female","Male")) +
  scale_x_continuous(limits=c(0,15.2), breaks=c(seq(0,15,2))) +
  labs(y="Gender", x="Predicted Composite Score (of 15)") +
  scale_color_manual("Site", labels=c("Choma","Mutasa","Nchelenge"),values = c(SteppedSequential5Steps[c(21,18,15)])) +
  theme_classic() +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=15),
        legend.title = element_text(size=14),
        legend.title.align=0.5,
        legend.text = element_text(size=12))
#
#


###### by trip # and trip time very similar to by all (trip #, trip time and clust time), where loc significant ###### 
###### but mainly due to slightly lower predicted value for Nchelenge (8.5 vs 10.2 and 10.7 for C and M) (out of 16 --> range 4-16) ###### 
###### similarly for all, Nchelenge has 11.0 vs 12.7 and 13.1 for c and M (out of 20 --> range 5-20) ###### 
###### by clust time (without trip count or time) is different because Choma higher (8.7 vs 7.0 and 7.1 for M and N) (out of 12 --> range 3-12) ###### 
###
###### by most (with trip count but not time since otherwise effect of trip # doubled) 
###### has Choma highest, Mutasa middle and Nchelenge lowest ( 10.8 vs 9.9 vs 9.3) (out of 16 --> range 4-16) ###### 


###             Choma     Mutasa      Nchelenege
### trip time   10.2       10.7          8.5    of 16
### clust time   8.7        7.0          7.1    of 12
### all         12.7       13.1          11.0   of 20
### most        10.8        9.9          9.3    of 16

### trip time    7.2        7.7          5.5    of 1-13
### clust time   6.7        5.0          5.1    of 1-10
### all          8.7        9.1          7.0    of 1-16
### most         7.8        6.9          6.3    of 1-13

###
##
##
##
########
#####
#####
#####
#####
#####