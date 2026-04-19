
# install.packages("ncdf4")

setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/R_Data/Livneh_precip/")

library(ncdf4)
library(lubridate)
library(tidyverse)

my.full.dataframe <- data.frame(Date = POSIXct(), lon = numeric(), lat = numeric(), variable_value = numeric())


for(i in 1915:2018){
  
  url <- paste0("https://psl.noaa.gov/thredds/dodsC/Datasets/livneh/metvars/altprecip/prec.", i, ".nc")
  nc <- nc_open(url)
  
  # print(nc)
  
  variable_array <- ncvar_get(nc, varid="prec")
  longitude <- ncvar_get(nc, varid="lon")
  latitude <- ncvar_get(nc, varid="lat")
  time <- ncvar_get(nc, varid="time")
  
  nc_close(nc)
  
  # Create a grid of all combinations of lon, lat, and time
  # This assumes 'variable_data' is a 3D array matching the dimensions
  indices <- expand.grid(lon = longitude, lat = latitude, time = time)
  
  # Combine the coordinate grid with the flattened (as.vector) variable data
  df <- data.frame(
    lon = indices$lon,
    lat = indices$lat,
    time = indices$time,
    variable_value = as.vector(variable_array)
  )
  
  # unique(df$lat) # 32.88
  # unique(df$lon) # 242.75
  
  df <- df[which(abs(df$lat - 32) < 1 & abs(df$lon - 242) < 1),]
  
  
  df$diff.lat <- abs(df$lat - 32.88)
  df$diff.lon <- abs(df$lon - 242.75)
  
  min.diff.lat <- unique(df$diff.lat)[order(unique(df$diff.lat), decreasing = F)][1:2]
  min.diff.lon <- unique(df$diff.lon)[order(unique(df$diff.lon), decreasing = F)][1:2]
  
  my.df <- df[which(df$diff.lat %in% min.diff.lat & df$diff.lon %in% min.diff.lon),]
  my.df$Date <- parse_date_time("1915-01-01", orders = "Ymd") + my.df$time * 60*60*24
  
  my.df <- my.df[,c("lon", "lat", "variable_value", "Date")]
  
  my.df <- my.df %>% group_by(Date) %>% summarise(across(everything(), ~mean(., na.rm = TRUE)))
  
  my.full.dataframe <- rbind(my.full.dataframe, my.df)
  
  print(paste0("Year ", i, " complete"))
  
}

saveRDS(my.full.dataframe, file = "2026-03-17_full_precip_df.rds")




