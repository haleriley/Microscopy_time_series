# quick personnel plot #
# RJH
# 2026-04-26

library(readxl)
library(tidyverse)
library(lubridate)


setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/")

personnel <- read_excel("R_Data/Microscopy_Carter/personnel.xlsx")

personnel$Year <- parse_date_time(personnel$Year, orders = "Y")
personnel$Year.end <- parse_date_time(as.numeric(year(personnel$Year))+1, orders = "Y")
personnel$Counter <- factor(personnel$Counter, levels = c("MH", "MC & MH", "MC", "MC & KS", "KS"))

ggplot(data = personnel) +
  geom_rect(aes(xmin = Year, xmax = Year.end, ymin = 0, ymax = 1, fill = Counter)) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  scale_fill_manual(values = scico(6, palette = "romaO", direction = -1)) +
  labs(fill = "", x = "Year") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), panel.grid.major.y = element_blank())


