# comparing microscopy-predicted O2bio to anomaly events (heatwaves, ENSO, etc)
# 2026-02-09
# RJH

# ---- library ----

library(tidyverse)
library(lubridate)
library(janitor)
library(vegan)
library(plotly)
library(readxl)
library(patchwork)
library(DESeq2)
library(akima)


setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/")

combo.full <- readRDS("2026-02-09_microscopy_o2bio_combo_full.rds")

enso.soi <- readRDS("R_Data/2026-02-09_enso_SOI.rds")
enso.sst <- readRDS("R_Data/2026-02-09_enso_SST.rds")

sio.temp <- readRDS("R_Data/2026-02-09_SIO_temp.rds")

# ---- smooth/average SIO temperature data ----

sio.temp.smooth <- sio.temp[1,]

combo <- merge(sio.temp, combo.full, by = "Date")

for(d in 1:nrow(combo)){
  
  my.date <- combo$Date[d]
  
  my.date1 <- my.date - 60*60*24
  my.date2 <- my.date + 60*60*47
  
  my.df <- sio.temp[which(sio.temp$Date >= my.date1 & sio.temp$Date <= my.date2),]
  
  my.df <- my.df %>% summarize_all(mean)
  
  my.df$Date <- my.date
  
  sio.temp.smooth <- rbind(sio.temp.smooth, my.df)
  
}

sio.temp.smooth <- sio.temp.smooth[-1,]

combo <- merge(sio.temp.smooth, combo.full, by = "Date")

# ---- compare surface water temp to predicted o2bio ----

ggplot(data = combo) +
  geom_point(aes(x = SURF_TEMP_C, y = O2bio.predicted.altered)) +
  geom_smooth(aes(x = SURF_TEMP_C, y = O2bio.predicted.altered), method = "lm", se = F, color = "red") +
  labs(x = expression("Surface Water Temp. (" * degree * "C)"), y = "Predicted [O2]bio (uM)") +
  ylim(c(-60,60)) +
  theme_bw() +
  annotate(geom = "text", x = 19, y = 50, label = "***", size = 8, color = "red")

summary(lm(combo$O2bio.predicted.altered~combo$SURF_TEMP_C))


# ---- compare to ENSO anomaly -----

## only starts in 1950/1951ish, so probably can't use this to compare across the full time series...
## maybe when we get the data back to 1992...

combo.full$Month <- month(combo.full$Date)
combo.full$Year <- year(combo.full$Date)

## SOI ##

combo.monthly <- combo.full %>% group_by(Year, Month) %>% summarize_all(mean)

colnames(enso.soi)[which(colnames(enso.soi) == "YEAR")] <- "Year"
colnames(enso.soi)[which(colnames(enso.soi) == "month")] <- "Month"
enso.soi$Month <- as.integer(factor(tools::toTitleCase(tolower(enso.soi$Month)), levels = month.abb))
enso.soi$Month[is.na(enso.soi$Month)] <- 5

combo.monthly <- merge(combo.monthly, enso.soi, by = c("Year", "Month"))

ggplot(data = combo.monthly) +
  geom_boxplot(aes(x = extreme..1, y = O2bio.predicted.altered)) +
  theme_bw()

t.test(x = combo.monthly$O2bio.predicted.altered[which(combo.monthly$extreme..1 == "yes")], y = combo.monthly$O2bio.predicted.altered[which(combo.monthly$extreme..1 == "no")])


## SST ##

combo.monthly <- combo.full %>% group_by(Year, Month) %>% summarize_all(mean)

colnames(enso.sst)[which(colnames(enso.sst) == "YR")] <- "Year"
colnames(enso.sst)[which(colnames(enso.sst) == "month")] <- "Month"
enso.sst$Month <- as.integer(factor(tools::toTitleCase(tolower(enso.sst$Month)), levels = month.abb))
enso.sst$Month[is.na(enso.sst$Month)] <- 5

combo.monthly <- merge(combo.monthly, enso.sst, by = c("Year", "Month"))

ggplot(data = combo.monthly) +
  geom_boxplot(aes(x = extreme..1, y = O2bio.predicted.altered)) +
  theme_bw()

t.test(x = combo.monthly$O2bio.predicted.altered[which(combo.monthly$extreme..1 == "yes")], y = combo.monthly$O2bio.predicted.altered[which(combo.monthly$extreme..1 == "no")])


# # ---- define SIO heatwaves (standard deviation) -- old methods ----
# 
# sio.temp.summary <- sio.temp %>% group_by(MONTH, DAY) %>% summarize(SURF_TEMP_C_MEAN = mean(SURF_TEMP_C, na.rm = T), SURF_TEMP_C_SD = sd(SURF_TEMP_C, na.rm = T))
# 
# sio.temp <- merge(sio.temp, sio.temp.summary, by = c("MONTH", "DAY"))
# 
# sio.temp$temp_high_1sd <- "no"
# sio.temp$temp_high_1sd[which(sio.temp$SURF_TEMP_C > (sio.temp$SURF_TEMP_C_MEAN + sio.temp$SURF_TEMP_C_SD))] <- "yes"
# 
# sio.temp$temp_high_2sd <- "no"
# sio.temp$temp_high_2sd[which(sio.temp$SURF_TEMP_C > (sio.temp$SURF_TEMP_C_MEAN + 2*sio.temp$SURF_TEMP_C_SD))] <- "yes"
# 
# sio.temp$temp_low_1sd <- "no"
# sio.temp$temp_low_1sd[which(sio.temp$SURF_TEMP_C < (sio.temp$SURF_TEMP_C_MEAN - sio.temp$SURF_TEMP_C_SD))] <- "yes"
# 
# sio.temp$temp_low_2sd <- "no"
# sio.temp$temp_low_2sd[which(sio.temp$SURF_TEMP_C < (sio.temp$SURF_TEMP_C_MEAN - 2*sio.temp$SURF_TEMP_C_SD))] <- "yes"
# 
# 
# sio.temp <- sio.temp[order(sio.temp$Date, decreasing = F),]
# 
# 
# sio.temp$heatwave_1sd <- "no"
# for(d in 2:(nrow(sio.temp)-1)){
#   
#   if(sio.temp$temp_high_1sd[d] == "yes" & sio.temp$temp_high_1sd[d-1] == "yes" & sio.temp$temp_high_1sd[d+1] == "yes"){
#    
#      sio.temp$heatwave_1sd[d] <- "yes"
#     sio.temp$heatwave_1sd[d-1] <- "yes"
#     sio.temp$heatwave_1sd[d+1] <- "yes"
#     
#   }
# }
# 
# sio.temp$heatwave_2sd <- "no"
# for(d in 2:(nrow(sio.temp)-1)){
#   
#   if(sio.temp$temp_high_2sd[d] == "yes" & sio.temp$temp_high_2sd[d-1] == "yes" & sio.temp$temp_high_2sd[d+1] == "yes"){
#     
#     sio.temp$heatwave_2sd[d] <- "yes"
#     sio.temp$heatwave_2sd[d-1] <- "yes"
#     sio.temp$heatwave_2sd[d+1] <- "yes"
#     
#   }
# }
# 
# sio.temp$coldwave_1sd <- "no"
# for(d in 2:(nrow(sio.temp)-1)){
#   
#   if(sio.temp$temp_low_1sd[d] == "yes" & sio.temp$temp_low_1sd[d-1] == "yes" & sio.temp$temp_low_1sd[d+1] == "yes"){
#     
#     sio.temp$coldwave_1sd[d] <- "yes"
#     sio.temp$coldwave_1sd[d-1] <- "yes"
#     sio.temp$coldwave_1sd[d+1] <- "yes"
#     
#   }
# }
# 
# sio.temp$coldwave_2sd <- "no"
# for(d in 2:(nrow(sio.temp)-1)){
#   
#   if(sio.temp$temp_low_2sd[d] == "yes" & sio.temp$temp_low_2sd[d-1] == "yes" & sio.temp$temp_low_2sd[d+1] == "yes"){
#     
#     sio.temp$coldwave_2sd[d] <- "yes"
#     sio.temp$coldwave_2sd[d-1] <- "yes"
#     sio.temp$coldwave_2sd[d+1] <- "yes"
#     
#   }
# }
# 
# 
# length(which(sio.temp$heatwave_1sd == "yes"))/nrow(sio.temp)
# length(which(sio.temp$heatwave_2sd == "yes"))/nrow(sio.temp)
# length(which(sio.temp$coldwave_1sd == "yes"))/nrow(sio.temp)
# length(which(sio.temp$coldwave_2sd == "yes"))/nrow(sio.temp)
# 
# 
# combo.sio.temp <- merge(combo.full, sio.temp, by = "Date")
# 
# 
# ggplot(data = combo.sio.temp) +
#   geom_boxplot(aes(x = heatwave_1sd, y = O2bio.predicted.altered)) +
#   theme_bw()
# 
# t.test(x = combo.sio.temp$O2bio.predicted.altered[which(combo.sio.temp$heatwave_1sd == "yes")], y = combo.sio.temp$O2bio.predicted.altered[which(combo.sio.temp$heatwave_1sd == "no")])
# 
# ggplot(data = combo.sio.temp) +
#   geom_boxplot(aes(x = heatwave_2sd, y = O2bio.predicted.altered)) +
#   theme_bw()
# 
# t.test(x = combo.sio.temp$O2bio.predicted.altered[which(combo.sio.temp$heatwave_2sd == "yes")], y = combo.sio.temp$O2bio.predicted.altered[which(combo.sio.temp$heatwave_2sd == "no")])
# 
# 
# ggplot(data = combo.sio.temp) +
#   geom_boxplot(aes(x = coldwave_1sd, y = O2bio.predicted.altered)) +
#   theme_bw()
# 
# t.test(x = combo.sio.temp$O2bio.predicted.altered[which(combo.sio.temp$coldwave_1sd == "yes")], y = combo.sio.temp$O2bio.predicted.altered[which(combo.sio.temp$coldwave_1sd == "no")])
# 
# ggplot(data = combo.sio.temp) +
#   geom_boxplot(aes(x = coldwave_2sd, y = O2bio.predicted.altered)) +
#   theme_bw()
# 
# t.test(x = combo.sio.temp$O2bio.predicted.altered[which(combo.sio.temp$coldwave_2sd == "yes")], y = combo.sio.temp$O2bio.predicted.altered[which(combo.sio.temp$coldwave_2sd == "no")])
# 

# ---- define SIO heatwaves (following Hobday et al. 2016) ----

## they recommend a 30-year time period for defining baseline climatologies
## I'm not sure if this is a minimum (i.e. at least 30 years) or an exact goal
## I am going to try a few ways of determining baseline climatologies

## Zhang et al (2005) recommends OOB bootstrapping but I don't know if this would be helpful here
## going to try and do a similar one 


sio.temp$cal.date <- parse_date_time(paste(sio.temp$MONTH, sio.temp$DAY, sep = "-"), orders = "md")

unique.cal.dates <- unique(sio.temp$cal.date)
unique.cal.dates <- unique.cal.dates[order(unique.cal.dates, decreasing = F)]

# sio.temp.baseline <- sio.temp[1:length(unique.cal.dates),]
sio.temp.baseline.bootstrap <- as.data.frame(unique.cal.dates)
colnames(sio.temp.baseline.bootstrap) <- "cal.date"
sio.temp.baseline.bootstrap <- sio.temp.baseline.bootstrap %>% arrange(cal.date)
  
## increase i eventually to make more comphrehensive bootstrapping
for(i in 1:100){
  
  my.years.boot <- sample(unique(sio.temp$YEAR), size = 1, replace = T)
  my.years.boot <- seq(my.years.boot, by = 1, length.out = 30)
  
  my.bootstrapped.df <- sio.temp[1,]
  my.bootstrapped.df <- my.bootstrapped.df[-1,]
  
  for(y in 1:length(my.years.boot)){
    
    my.bootstrapped.df <- rbind(my.bootstrapped.df, sio.temp[which(sio.temp$YEAR == my.years.boot[y]),])
    
  }
  
  
  my.baselines <- as.data.frame(matrix(nrow = length(unique.cal.dates), ncol = 2))
  colnames(my.baselines) <- c("cal.date", "boot.base.temp")
  my.baselines$cal.date <- unique.cal.dates
  
  for(d in 1:length(unique.cal.dates)){
    
    my.date <- unique.cal.dates[d]
    
    my.date.5before <- my.date - (60 * 60 * 24 *5)
    my.date.5after <- my.date + (60 * 60 * 24 * 5)
    
    year(my.date.5before) <- 0000
    year(my.date.5after) <- 0000
    
    if(my.date.5before > my.date.5after){
      
      my.df <- my.bootstrapped.df[which(my.bootstrapped.df$cal.date >= my.date.5before | my.bootstrapped.df$cal.date <= my.date.5after),]
      
      my.baseline <- my.df %>% summarize(temp.mean.boot = mean(SURF_TEMP_C, na.rm = T))
      
      my.baselines$boot.base.temp[which(my.baselines$cal.date == unique.cal.dates[d])] <- my.baseline[1,1]
      
    }
    
    if(my.date.5before < my.date.5after){
      
      my.df <- my.bootstrapped.df[which(my.bootstrapped.df$cal.date >= my.date.5before & my.bootstrapped.df$cal.date <= my.date.5after),]
      
      my.baseline <- my.df %>% summarize(temp.mean.boot = mean(SURF_TEMP_C, na.rm = T))
      
      my.baselines$boot.base.temp[which(my.baselines$cal.date == unique.cal.dates[d])] <- my.baseline[1,1]
      
    }
    
  }
  
  
  sio.temp.baseline.bootstrap <- merge(sio.temp.baseline.bootstrap, my.baselines, by = "cal.date", all = T)
  colnames(sio.temp.baseline.bootstrap)[i+1] <- paste0("boot.base.temp.", i)
  
  print(paste0("Bootstrap ", i, " of 100"))
  
}

sio.temp.baseline.bootstrap.long <- sio.temp.baseline.bootstrap %>% pivot_longer(cols = colnames(sio.temp.baseline.bootstrap[-1]), names_to = "iteration", values_to = "boot.base.temp")
sio.temp.baseline.bootstrap.summary <- sio.temp.baseline.bootstrap.long %>% group_by(cal.date) %>% dplyr::summarize(boot.base.temp = mean(boot.base.temp, na.rm = T))

# 
# for(d in 1:length(unique.cal.dates)){
#   
#   my.date <- unique.cal.dates[d]
#   
#   my.date.5before <- my.date - (60 * 60 * 24 *5)
#   my.date.5after <- my.date + (60 * 60 * 24 * 5)
#   
#   year(my.date.5before) <- 0000
#   year(my.date.5after) <- 0000
#   
#   if(my.date.5before > my.date.5after){
#     
#     my.df <- sio.temp[which(sio.temp$cal.date >= my.date.5before | sio.temp$cal.date <= my.date.5after),]
#     
#     sio.temp.baseline[d,] <- my.df %>% summarize_all(mean, na.rm = T)
#     sio.temp.baseline$cal.date[d] <- my.date
#     
#   }
#   
#   if(my.date.5before < my.date.5after){
#     
#     my.df <- sio.temp[which(sio.temp$cal.date >= my.date.5before & sio.temp$cal.date <= my.date.5after),]
#     
#     sio.temp.baseline[d,] <- my.df %>% summarize_all(mean, na.rm = T)
#     sio.temp.baseline$cal.date[d] <- my.date
#     
#   }
#   
# }


ggplot(data = sio.temp.baseline.bootstrap.summary) +
  geom_line(aes(x = cal.date, y = boot.base.temp)) +
  theme_bw()

sio.temp.baseline <- sio.temp.baseline.bootstrap.summary

colnames(sio.temp.baseline)[which(colnames(sio.temp.baseline) == "boot.base.temp")] <- "SURF_TEMP_C_BASELINE"

sio.temp <- merge(sio.temp, sio.temp.baseline[,c("cal.date", "SURF_TEMP_C_BASELINE")], by = "cal.date")

sio.temp$SURF_TEMP_ANOM <- sio.temp$SURF_TEMP_C - sio.temp$SURF_TEMP_C_BASELINE

ggplot(data = sio.temp) +
  geom_line(aes(x = Date, y = SURF_TEMP_ANOM)) +
  theme_bw()


sio.temp$percentile.90.hot <- "no"
sio.temp$percentile.90.hot[which(sio.temp$SURF_TEMP_ANOM >= quantile(sio.temp$SURF_TEMP_ANOM[which(sio.temp$SURF_TEMP_ANOM > 0)], 0.9, na.rm = T))] <- "yes"

sio.temp$percentile.90.cold <- "no"
sio.temp$percentile.90.cold[which(sio.temp$SURF_TEMP_ANOM <= -1*quantile(abs(sio.temp$SURF_TEMP_ANOM[which(sio.temp$SURF_TEMP_ANOM < 0)]), 0.9, na.rm = T))] <- "yes"

## toggle this if only looking for heatwave/coldwaves outside of range of climatological baseline (highest ecological stress)
# baseline.min <- min(sio.temp$SURF_TEMP_C_BASELINE)
# baseline.max <- max(sio.temp$SURF_TEMP_C_BASELINE)
# 
# sio.temp$percentile.90.hot[which(sio.temp$SURF_TEMP_C < baseline.max)] <- "no"
# sio.temp$percentile.90.cold[which(sio.temp$SURF_TEMP_C > baseline.min)] <- "no"


sio.temp$heatwave.perc.90 <- "no"
sio.temp$coldwave.perc.90 <- "no"


sio.temp <- sio.temp[order(sio.temp$Date, decreasing = F),]


## define heatwaves as periods of 5 days or more of 90th-percentile anomaly
for(d in 3:(nrow(sio.temp)-2)){

  if(length(which(sio.temp$percentile.90.hot[c(d-2,d-1,d,d+1,d+2)] == "yes")) == 5){

    sio.temp$heatwave.perc.90[c(d-2,d-1,d,d+1,d+2)] <- "yes"

  }
  
  if(length(which(sio.temp$percentile.90.cold[c(d-2,d-1,d,d+1,d+2)] == "yes")) == 5){
    
    sio.temp$coldwave.perc.90[c(d-2,d-1,d,d+1,d+2)] <- "yes"
    
  }
  
}

## combine heatwaves with only one or two "yes heatwave" days betweeen them
for(d in 3:(nrow(sio.temp)-2)){
  
  if(length(which(sio.temp$heatwave.perc.90[c(d-2,d-1,d,d+1,d+2)] == "yes")) >= 3){
    
    sio.temp$heatwave.perc.90[d] <- "yes"
    
  }
  
  if(length(which(sio.temp$coldwave.perc.90[c(d-2,d-1,d,d+1,d+2)] == "yes")) >= 3){
    
    sio.temp$coldwave.perc.90[d] <- "yes"
    
  }
  
}


# ---- plot SIO water temp with heat and coldwaves ----

# p <- ggplot(sio.temp) +
#   geom_line(aes(x = Date, y = SURF_TEMP_C)) +
#   geom_line(aes(x = Date, y = SURF_TEMP_C_BASELINE), color = "gold") +
#   theme_bw()
# ggplotly(p)


# ggplot() +
#   geom_line(data = sio.temp, aes(x = Date, y = (SURF_TEMP_C - SURF_TEMP_C_MEAN)/SURF_TEMP_C_SD), alpha = 0.4) +
#   geom_point(data = sio.temp[which(sio.temp$heatwave_2sd == "yes"),], aes(x = Date, y = (SURF_TEMP_C - SURF_TEMP_C_MEAN)/SURF_TEMP_C_SD), color = "red") +
#   geom_point(data = sio.temp[which(sio.temp$coldwave_2sd == "yes"),], aes(x = Date, y = (SURF_TEMP_C - SURF_TEMP_C_MEAN)/SURF_TEMP_C_SD), color = "blue") +
#   # geom_point(data = sio.temp[which(sio.temp$heatwave_1sd == "yes"),], aes(x = Date, y = (SURF_TEMP_C - SURF_TEMP_C_MEAN)/SURF_TEMP_C_SD), color = "red", alpha = 0.4) +
#   # geom_point(data = sio.temp[which(sio.temp$coldwave_1sd == "yes"),], aes(x = Date, y = (SURF_TEMP_C - SURF_TEMP_C_MEAN)/SURF_TEMP_C_SD), color = "blue", alpha = 0.4) +
#   theme_bw()


ggplot() +
  geom_line(data = sio.temp, aes(x = Date, y = SURF_TEMP_C), alpha = 0.4) +
  geom_smooth(data = sio.temp, aes(x = Date, y = SURF_TEMP_C), color = "purple", method = "lm", se = F) +
  geom_point(data = sio.temp[which(sio.temp$heatwave.perc.90 == "yes"),], aes(x = Date, y = SURF_TEMP_C, color = "Heatwave"), alpha = 1) +
  geom_point(data = sio.temp[which(sio.temp$coldwave.perc.90 == "yes"),], aes(x = Date, y = SURF_TEMP_C, color = "Coldwave"), alpha = 1) +
  labs(y = expression("Surface Water Temp. (" * degree * "C)")) +
  theme_bw() +
  scale_color_manual(name = "", values = c("Heatwave" = "red", "Coldwave" = "blue"), labels = rev(factor(c("Heatwave", "Coldwave"), levels = c("Heatwave", "Coldwave"))))

hist(sio.temp$MONTH[which(sio.temp$heatwave.perc.90 == "yes")], 
     main = "Heatwave (90th percentile) Months Histogram",
     xlab = "Month")

hist(sio.temp$MONTH[which(sio.temp$coldwave.perc.90 == "yes")], 
     main = "Coldwave (90th percentile) Months Histogram",
     xlab = "Month")


ggplot(data = sio.temp) +
  geom_violin(aes(x = MONTH, y = SURF_TEMP_C, group = MONTH)) +
  theme_bw()


length(which(sio.temp$heatwave.perc.90 == "yes"))/nrow(sio.temp)
length(which(sio.temp$coldwave.perc.90 == "yes"))/nrow(sio.temp)

combo.sio.temp <- merge(combo.full, sio.temp, by = "Date", all = T)
combo.sio.temp.overlap <- merge(combo.full, sio.temp, by = "Date", all = F)

o2bio.interp <- as.data.frame(approx(x = combo.sio.temp.overlap$Date, y = combo.sio.temp.overlap$O2bio.predicted.altered, xout = combo.sio.temp$Date))
colnames(o2bio.interp) <- c("Date", "O2bio.predicted.altered.interp")
               

combo.sio.temp <- merge(x = combo.sio.temp, y = o2bio.interp, by = "Date", all.x = T)


h <- ggplot(data = combo.sio.temp) +
  geom_violin(aes(x = heatwave.perc.90, y = O2bio.predicted.altered.interp, fill = heatwave.perc.90)) +
  labs(x = "Heatwave", y = "Predicted [O2]bio") +
  scale_x_discrete(labels = c("No", "Yes")) +
  scale_fill_manual(values = c("grey", "red")) +
  theme_bw() +
  theme(legend.position = "none") +
  annotate("text", y = 40, x = "yes", label = "***", size = 8)

t.test(x = combo.sio.temp$O2bio.predicted.altered.interp[which(combo.sio.temp$heatwave.perc.90 == "yes")], y = combo.sio.temp$O2bio.predicted.altered.interp[which(combo.sio.temp$heatwave.perc.90 == "no")])


c <- ggplot(data = combo.sio.temp) +
  geom_violin(aes(x = coldwave.perc.90, y = O2bio.predicted.altered.interp, fill = coldwave.perc.90)) +
  labs(x = "Coldwave", y = "Predicted [O2]bio") +
  scale_x_discrete(labels = c("No", "Yes")) +
  scale_fill_manual(values = c("grey", "blue")) +
  theme_bw() +
  theme(legend.position = "none") +
  annotate("text", y = 40, x = "yes", label = "***", size = 8)


t.test(x = combo.sio.temp$O2bio.predicted.altered.interp[which(combo.sio.temp$coldwave.perc.90 == "yes")], y = combo.sio.temp$O2bio.predicted.altered.interp[which(combo.sio.temp$coldwave.perc.90 == "no")])

h+c

p <- ggplot(data = combo.sio.temp) +
  geom_line(aes(x = Date, y = O2bio.predicted.altered.interp), alpha = 0.4) +
  geom_point(data = combo.sio.temp[which(combo.sio.temp$heatwave.perc.90 == "yes"),], aes(x = Date, y = O2bio.predicted.altered.interp), color = "red") +
  geom_point(data = combo.sio.temp[which(combo.sio.temp$coldwave.perc.90 == "yes"),], aes(x = Date, y = O2bio.predicted.altered.interp), color = "blue") +
  theme_bw()
ggplotly(p)


saveRDS(combo.sio.temp, file = "2026-04-14_combo_sio_temp.rds")


# ---- investigate correlation between o2bio change/slope and water temp (NONE) ----

combo.sio.temp <- combo.sio.temp[which(is.na(combo.sio.temp$O2bio.predicted.altered.interp) == F),]

my.lm.cor.df <- as.data.frame(matrix(nrow = nrow(combo.sio.temp), ncol = 4))
colnames(my.lm.cor.df) <- c("slope", "p", "R2", "RMSE")

for(d in 3:(nrow(combo.sio.temp)-2)){
  
  my.lm <- lm(combo.sio.temp$O2bio.predicted.altered.interp[(d-2):(d+2)]~combo.sio.temp$Date[(d-2):(d+2)])
  
  my.lm.cor.df$slope[d] <- summary(my.lm)$coefficients[2,1]
  my.lm.cor.df$p[d] <- summary(my.lm)$coefficients[2,4]
  my.lm.cor.df$R2[d] <- summary(my.lm)$adj.r.squared
  my.lm.cor.df$RMSE[d] <- sqrt(mean((summary(my.lm)$residuals)^2))
  
}

combo.sio.temp <- cbind(combo.sio.temp, my.lm.cor.df)


ggplot(data = combo.sio.temp) +
  geom_boxplot(aes(x = heatwave.perc.90, y = slope)) +
  theme_bw()

t.test(x = combo.sio.temp$slope[which(combo.sio.temp$heatwave.perc.90 == "yes")], y = combo.sio.temp$slope[which(combo.sio.temp$heatwave.perc.90 == "no")])

ggplot(data = combo.sio.temp) +
  geom_boxplot(aes(x = coldwave.perc.90, y = slope)) +
  theme_bw()

t.test(x = combo.sio.temp$slope[which(combo.sio.temp$coldwave.perc.90 == "yes")], y = combo.sio.temp$slope[which(combo.sio.temp$coldwave.perc.90 == "no")])

summary(lm(combo.sio.temp$slope~combo.sio.temp$SURF_TEMP_C))
ggplot(data = combo.sio.temp) +
  geom_point(aes(x = SURF_TEMP_C, y = slope)) +
  geom_smooth(aes(x = SURF_TEMP_C, y = slope), method = "lm", se = F) +
  theme_bw()













# ----






















