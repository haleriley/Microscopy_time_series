# extreme events data curation
# RJH
# 2026-01-27

# ---- library ----

library(lubridate)
library(tidyverse)


# ---- read in raw data ----

setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/R_Data/")

enso.sst <- read.table(file = "2026-01-27_ENSO_ONI3.4_SST_anomaly.txt", header = T)
enso.soi <- read.table(file = "2026-01-27_ENSO_SOI.txt", header = T, skip = 87, nrows = 75)

storms <- read.csv("storm_data_search_results.csv")

sio.temp <- read.csv("SIO_shore_stations_SST_sal/LaJolla_TEMP_1916-202509.csv", header = T, skip = 46)

# ---- clean data ----

## ENSO SST ##

enso.sst$month <- rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), length.out = nrow(enso.sst))
enso.sst$Date <- parse_date_time(paste(enso.sst$month, enso.sst$YR, sep = "-"), orders = "bY")
enso.sst$extreme..0.5 <- "no"
enso.sst$extreme..1 <- "no"

for(m in 3:(nrow(enso.sst)-2)){
  
  my.df <- enso.sst[(m-2):(m+2),]
  
  if(length(which(abs(my.df$ANOM) >= 0.5))){
    enso.sst$extreme..0.5[m] <- "yes"
  }
  
  if(length(which(abs(my.df$ANOM) >= 1))){
    enso.sst$extreme..1[m] <- "yes"
  }
}

## ENSO SOI ###

enso.soi <- enso.soi %>% pivot_longer(cols = colnames(enso.soi)[2:13], names_to = "Month", values_to = "standardized.anomaly")
enso.soi$Date <- parse_date_time(paste(enso.soi$YEAR, enso.soi$Month, sep = "-"), orders = "Yb")

enso.soi$extreme..0.5 <- "no"
enso.soi$extreme..1 <- "no"

for(m in 3:(nrow(enso.soi)-2)){
  
  my.df <- enso.soi[(m-2):(m+2),]
  
  if(length(which(abs(my.df$standardized.anomaly) >= 0.5))){
    enso.soi$extreme..0.5[m] <- "yes"
  }
  
  if(length(which(abs(my.df$standardized.anomaly) >= 1))){
    enso.soi$extreme..1[m] <- "yes"
  }
}


## SIO TEMP ##

sio.temp$Date <- parse_date_time(paste(sio.temp$YEAR, sio.temp$MONTH, sio.temp$DAY, sep = "-"), orders = "Ymd")



# ---- visualize data ----

ggplot(data = enso.sst) +
  geom_hline(yintercept = 0.5, color = "red", linewidth = 1, alpha = 0.5) +
  geom_hline(yintercept = -0.5, color = "blue", linewidth = 1, alpha = 0.5) +
  geom_hline(yintercept = 1, color = "darkred", linewidth = 1, alpha = 0.5) +
  geom_hline(yintercept = -1, color = "darkblue", linewidth = 1, alpha = 0.5) +
  geom_line(aes(x = Date, y = ANOM)) +
  labs(y = "SST Anomaly") +
  theme_bw()


ggplot(data = enso.soi) +
  geom_hline(yintercept = 0.5, color = "red", linewidth = 1, alpha = 0.5) +
  geom_hline(yintercept = -0.5, color = "blue", linewidth = 1, alpha = 0.5) +
  geom_hline(yintercept = 1, color = "darkred", linewidth = 1, alpha = 0.5) +
  geom_hline(yintercept = -1, color = "darkblue", linewidth = 1, alpha = 0.5) +
  geom_line(aes(x = Date, y = standardized.anomaly)) +
  labs(y = "Standardized SOI Anomaly") +
  theme_bw()


ggplot(data = sio.temp) +
  geom_line(aes(x = Date, y = SURF_TEMP_C)) +
  theme_bw()

# ---- save curated data ----

saveRDS(enso.soi, file = "2026-02-09_enso_SOI.rds")
saveRDS(enso.sst, file = "2026-02-09_enso_SST.rds")
saveRDS(sio.temp, file = "2026-02-09_SIO_temp.rds")







