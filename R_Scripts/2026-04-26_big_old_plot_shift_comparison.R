# analysis of parameters for big ole plot and comparisons between shifts
# RJH
# 2026-04-26




# ---- library ----

library(tidyverse)
library(lubridate)
library(janitor)
library(vegan)
library(plotly)
library(readxl)
library(patchwork)
library(DESeq2)
library(gplots)
library(viridis)
library(ggthemes)
library(khroma)
library(ggsignif)
library(multcompView)
library(corrplot)
library(Hmisc)
library(betapart)



# ---- read in data ----

setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/R_Data/")

# combo.full <- readRDS("2026-02-09_microscopy_o2bio_combo_full.rds")
# combo.full <- readRDS("2026-03-23_microscopy_o2bio_combo_full.rds")
combo.full <- readRDS("2026-04-23_microscopy_o2bio_combo_full.rds")
# combo.sio.temp <- readRDS("2026-02-14_combo_sio_temp.rds") ## heatwaves only outside climatological baseline
combo.sio.temp <- readRDS("2026-02-14_combo_sio_temp.rds") ## all heatwaves

model.var.imp <- readRDS("2026-02-09_microscopy_model_var_imp.rds")

# predictors <- readRDS("2026-02-09_microscopy_model_predictors.rds")
# predictors <- readRDS("2026-03-23_phyto_rf_predictors.rds")
predictors <- readRDS("2026-03-27_phyto_rf_predictors.rds")
# predictors <- predictors[which(substr(predictors, 1, 5) != "other")]

unique.dinos <- readRDS("2026-04-23_unique_dinos.rds")
unique.diatoms <- readRDS("2026-04-23_unique_diatomss.rds")

model.var.imp <- readRDS("2026-02-09_microscopy_model_var_imp.rds")

o2bio.daily <- readRDS("../../O2-Ar_time_series/R_Data/2025-10-07_o2bio_est_df_daily.rds")


setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/")



phytos.combo <- readRDS("R_Data/2026-04-26_phytos_combo.rds")
phyto.mat <- readRDS("R_Data/2026-04-26_phyto_mat.rds")
combo.sio.temp <- readRDS("R_Data/2026-04-26_combo_sio_temp.rds")
combo.full <- readRDS("R_Data/2026-04-26_combo_full.rds")



bright <- color("bright")
my.colorblind.colors <- bright(7)
my.colors.blob.2022 <- NA

# ---- look at three shifts/four periods throughout time series - shifting taxa, NMDS Dim1, and 02bio values ----

## copied from above to remind:

heat.cols <- viridis::viridis(100, direction = 1)

phyto.mat.log10 <- log10(phyto.mat)
phyto.mat.log10[which(phyto.mat.log10 < 0)] <- 0
phyto.mat.log10 <- t(phyto.mat.log10)

summary(rowSums(phyto.mat.log10))
phyto.mat.log10 <- phyto.mat.log10[which(rowSums(phyto.mat.log10) >= 200),]
my.abund.taxa <- rownames(phyto.mat.log10)

my.colCol <- rep(my.colorblind.colors[4], times = ncol(phyto.mat.log10))
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) <= 2013)] <- my.colorblind.colors[3]
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) %in% c(2014, 2015))] <- my.colorblind.colors[2]
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) >= 2022)] <- my.colorblind.colors[6]


heatmap.2(phyto.mat.log10,
          trace = 'none',
          dendrogram = "row",
          scale = NULL,
          col = heat.cols,
          cexRow = 0.4,
          cexCol = 0.5,
          key = FALSE,
          margins = c(3, 8),
          srtCol = 90,
          lhei = c(1,100),
          offsetCol = 0,
          # adjCol = 0.5,
          Colv = F,
          colCol = my.colCol
          # key = TRUE,
          # keysize = 2,
          # density.info = "density"
)

# heatmap.2(phyto.mat.log10,
#           trace = 'none',
#           dendrogram = "row",
#           scale = NULL,
#           col = heat.cols,
#           cexRow = 0.4,
#           cexCol = 0.5,
#           key = FALSE,
#           margins = c(3, 8),
#           srtCol = 90,
#           lhei = c(1,100),
#           offsetCol = 0,
#           adjCol = 0.5,
#           Colv = T,
#           colCol = my.colCol
#           # key = TRUE,
#           # keysize = 2,
#           # density.info = "density"
# )


phyto.df <- as.data.frame(phyto.mat)
phyto.df.long <- phyto.df %>% pivot_longer(cols = colnames(phyto.df), names_to = "Taxon", values_to = "Count")

my.richness <- specnumber(phyto.df)
my.shannon.div <- diversity(phyto.df, index = "shannon")

phyto.df$richness <- my.richness
phyto.df$shannon.div <- my.shannon.div

phyto.df$Date <- parse_date_time(rownames(phyto.df), orders = "Ymd")

ggplot(phyto.df) +
  geom_line(aes(x = Date, y = richness)) +
  geom_smooth(aes(x = Date, y = richness), method = "loess") +
  labs(x = "Date", y = "Phytoplankton Richness") +
  theme_bw() +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") 


ggplot(phyto.df) +
  geom_line(aes(x = Date, y = shannon.div)) +
  geom_smooth(aes(x = Date, y = shannon.div), method = "loess") +
  labs(x = "Date", y = "Shannon Diversity") +
  theme_bw() +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") 


phyto.df$Blob.status <- "Blob"
phyto.df$Blob.status[which(year(phyto.df$Date) <= 2013)] <- "PreBlob"
phyto.df$Blob.status[which(year(phyto.df$Date) >= 2016)] <- "PostBlob"
phyto.df$Blob.status[which(year(phyto.df$Date) >= 2022)] <- "Post2022Shift"
phyto.df$Blob.status <- factor(phyto.df$Blob.status, levels = c("PreBlob", "Blob", "PostBlob", "Post2022Shift"))


## BEUTI ##

beuti.df <- read.csv("R_Data/BEUTI_daily.csv")
beuti.df$Date <- parse_date_time(paste(beuti.df$year, beuti.df$month, beuti.df$day, dep = "-"), orders = "Ymd")
beuti.df <- beuti.df[,which(colnames(beuti.df) %in% c("Date", "X33N"))]
colnames(beuti.df)[1] <- "BEUTI.33N"

phyto.df <- merge(phyto.df, beuti.df, by = "Date", all.x = T, all.y = F)

my.aov <- aov(phyto.df$BEUTI.33N~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$BEUTI.33N, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.beuti.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = BEUTI.33N), fill = "grey30") +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "BEUTI Index \n 33N \n (mmol/m/sec)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = log10(y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


## Richness ##
my.aov <- aov(phyto.df$richness~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$richness, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.richness.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = richness), fill = my.colorblind.colors[6]) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Phytoplankton \n Richness") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = log10(y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


## Shannon Diversity ##
my.aov <- aov(phyto.df$shannon.div~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$shannon.div, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.shannon.div.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = shannon.div), fill = my.colorblind.colors[1]) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Phytoplankton \n Diversity") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


## O2bio ##

phyto.df.o2biodates <- merge(phyto.df, o2bio.daily[,which(colnames(o2bio.daily) %in% c("Date", "O2bio.estimated"))], all = T)

phyto.df.o2biodates$Blob.status <- "Blob"
phyto.df.o2biodates$Blob.status[which(year(phyto.df.o2biodates$Date) <= 2013)] <- "PreBlob"
phyto.df.o2biodates$Blob.status[which(year(phyto.df.o2biodates$Date) >= 2016)] <- "PostBlob"
phyto.df.o2biodates$Blob.status[which(year(phyto.df.o2biodates$Date) >= 2022)] <- "Post2022Shift"
phyto.df.o2biodates$Blob.status <- factor(phyto.df.o2biodates$Blob.status, levels = c("PreBlob", "Blob", "PostBlob", "Post2022Shift"))


my.aov <- aov(phyto.df.o2biodates$O2bio.estimated~phyto.df.o2biodates$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df.o2biodates$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df.o2biodates$O2bio.estimated, phyto.df.o2biodates$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.o2bio.box <- ggplot(phyto.df.o2biodates) +
  geom_boxplot(aes(x = Blob.status, y = O2bio.estimated), fill = "black") +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "[O2]bio (uM)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))



## water temperature ##

phyto.df <- merge(phyto.df, combo.sio.temp[,which(colnames(combo.sio.temp) %in% c("Date", "SURF_TEMP_C"))])

my.aov <- aov(phyto.df$SURF_TEMP_C~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$SURF_TEMP_C, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.temp.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = SURF_TEMP_C), fill = my.colorblind.colors[2]) +
  # scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Surface Water \n Temperature (C)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = y_pos, label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

# ggplot(phyto.df) +
#   geom_line(aes(x = Date, y = SURF_TEMP_C)) +
#   geom_smooth(aes(x = Date, y = SURF_TEMP_C), method = "loess") +
#   labs(x = "Date", y = "Surface Water Temperature (C)") +
#   theme_bw() +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y")  +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


## water temperature anomaly ##

phyto.df <- merge(phyto.df, combo.sio.temp[,which(colnames(combo.sio.temp) %in% c("Date", "SURF_TEMP_ANOM"))])

my.aov <- aov(phyto.df$SURF_TEMP_ANOM~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$SURF_TEMP_ANOM, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.tempanom.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = SURF_TEMP_ANOM), fill = my.colorblind.colors[4]) +
  # scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Surface Water \n Temperature \n Anomaly (C)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos*2), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

# ggplot(phyto.df) +
#   geom_line(aes(x = Date, y = SURF_TEMP_ANOM)) +
#   geom_smooth(aes(x = Date, y = SURF_TEMP_ANOM), method = "loess") +
#   labs(x = "Date", y = "Surface Water Temperature Anomaly (C)") +
#   theme_bw() +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y")  +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


## phytoplankton total counts ##

phyto.df$sum.phyto.counts <- rowSums(phyto.df[,which(colnames(phyto.df) %in% predictors)])

my.aov <- aov(log10(phyto.df$sum.phyto.counts)~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(log10(phyto.df$sum.phyto.counts), phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.phytocount.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = log10(sum.phyto.counts)), fill = my.colorblind.colors[5]) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "log10 Total \n Phytoplankton \n Count") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

# ggplot(phyto.df) +
#   geom_line(aes(x = Date, y = log10(sum.phyto.counts))) +
#   geom_smooth(aes(x = Date, y = log10(sum.phyto.counts)), method = "loess") +
#   labs(x = "Date", y = "log10(Total Phytoplankton) (Counts)") +
#   theme_bw() +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y")  +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


## chlorophyll ##

sccoos <- read.csv("https://erddap.sccoos.org/erddap/tabledap/HABs-ScrippsPier.csv?time%2CChl_Volume_Filtered%2CChl1%2CChl2%2CAvg_Chloro%2CPhosphate%2CSilicate%2CNitrite%2CNitrite_Nitrate%2CAmmonium%2CNitrate%2CVolume_Settled_for_Counting&time%3E=2008-06-30T15%3A00%3A00Z&time%3C=2026-04-06T18%3A16%3A00Z")

sccoos <- sccoos[-1,]
colnames(sccoos)[1] <- c("Date")
sccoos$Date <- parse_date_time(sccoos$Date, orders= "Ymd HMS")
sccoos$Date <- parse_date_time(paste(year(sccoos$Date), month(sccoos$Date), day(sccoos$Date), sep = "-"), orders = "Ymd")
sccoos <- sccoos %>% mutate(across(c(colnames(sccoos)[c(2:ncol(sccoos))]), as.numeric))
sccoos <- sccoos %>% group_by(Date) %>% summarize_all(mean)
phyto.df <- merge(phyto.df, sccoos, by = "Date", all.x = T, all.y = F)

my.aov <- aov(log10(phyto.df$Avg_Chloro)~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(log10(phyto.df$Avg_Chloro), phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.chla.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = log10(Avg_Chloro)), fill = my.colorblind.colors[3]) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "log10 \n Chlorophyll \n (ug/L)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

# ggplot(phyto.df) +
#   geom_line(aes(x = Date, y = log10(Avg_Chloro))) +
#   geom_smooth(aes(x = Date, y = log10(Avg_Chloro)), method = "loess") +
#   labs(x = "Date", y = "log10(Chlorophyll) (ug/L)") +
#   theme_bw() +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y")  +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


## nitrate ##

my.aov <- aov(phyto.df$Nitrate~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$Nitrate, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.nitrate.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = Nitrate), fill = "white") +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Nitrate (uM)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

# ggplot(phyto.df) +
#   geom_line(aes(x = Date, y = Nitrate)) +
#   geom_smooth(aes(x = Date, y = Nitrate), method = "loess") +
#   labs(x = "Date", y = "Nitrate (uM)") +
#   theme_bw() +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y")  +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


## phosphate ##

my.aov <- aov(phyto.df$Phosphate~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$Phosphate, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.phosphate.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = Phosphate), fill = "white") +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Phosphate (uM)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

# ggplot(phyto.df) +
#   geom_line(aes(x = Date, y = Phosphate)) +
#   geom_smooth(aes(x = Date, y = Phosphate), method = "loess") +
#   labs(x = "Date", y = "Phosphate (uM)") +
#   theme_bw() +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y")  +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


## silicate ##

my.aov <- aov(phyto.df$Silicate~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$Silicate, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.silicate.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = Silicate), fill = "white") +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Silicate (uM)") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = (y_pos), label = Letter),
    inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

# ggplot(phyto.df) +
#   geom_line(aes(x = Date, y = Silicate)) +
#   geom_smooth(aes(x = Date, y = Silicate), method = "loess") +
#   labs(x = "Date", y = "Silicate (uM)") +
#   theme_bw() +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y")  +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


## Predicted O2bio ##

phyto.df <- merge(phyto.df, combo.full[,which(colnames(combo.full) %in% c("Date", "O2bio.predicted.altered"))], by = "Date")

my.aov <- aov(phyto.df$O2bio.predicted.altered~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$O2bio.predicted.altered, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


plot.o2bio.predalt.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = O2bio.predicted.altered), fill = "black") +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Predicted \n [O2]bio (uM)") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = (y_pos), label = Letter),
    inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))





## comparisons

ggplot(phyto.df) +
  geom_point(aes(x = SURF_TEMP_ANOM, y = log10(Avg_Chloro))) +
  geom_smooth(aes(x = SURF_TEMP_ANOM, y = log10(Avg_Chloro)), method = "lm") +
  labs(x = "Surface Water Temperature Anomaly (C)", y = "log10(Chlorophyll) (ug/L)") +
  theme_bw()
summary(lm(log10(phyto.df$Avg_Chloro)~phyto.df$SURF_TEMP_ANOM))





# ---- measure volatility ----



# sccoos.FCM <- readRDS("../SCCOOS_microbial_time_series/R_Data/2026-03-17_sccoos_com_df_fcm_cleaned_rarefied_2000.rds")
# sccoos.FCM <- sccoos.FCM[,which(colnames(sccoos.FCM) %in% c("Date", "total_ml_vol_corrected_AF", "total_ml_vol_corrected_SG"))]
# phyto.df <- merge(phyto.df, sccoos.FCM, by = "Date", all.x = T, all.y = F)

# phyto.df <- merge(phyto.df, combo.full[which(colnames(combo.full) %in% c("Date", "O2bio.estimated"))], all.x = T, all.y = F)

phyto.df$sum.phyto.counts <- rowSums(phyto.df[,which(colnames(phyto.df) %in% predictors)])



phyto.df$volatility.O2bio.predicted.altered <- NA
phyto.df$volatility.chlorophyll <- NA
phyto.df$volatility.fcm.AF <- NA
phyto.df$volatility.log10chlorophyll <- NA
phyto.df$volatility.sum.count <- NA
phyto.df$volatility.log10.sum.count <- NA
phyto.df$volatility.temperature <- NA
phyto.df$volatility.temperature.anomaly <- NA
phyto.df$volatility.richness <- NA
phyto.df$volatility.shannon.div <- NA
phyto.df$volatility.nitrate <- NA
phyto.df$volatility.phosphate <- NA
phyto.df$volatility.silicate <- NA
phyto.df$volatility.beuti <- NA



for(d in 1:nrow(phyto.df)){
  
  my.date <- phyto.df$Date[d]
  
  my.df <- phyto.df[which(phyto.df$Date >= (my.date - 60*60*24*15) & phyto.df$Date <= (my.date + 60*60*24*15)),]
  
  phyto.df$volatility.O2bio.predicted.altered[d] <- sd(my.df$O2bio.predicted.altered, na.rm = T)
  phyto.df$volatility.chlorophyll[d] <- sd(my.df$Avg_Chloro, na.rm = T)
  phyto.df$volatility.log10chlorophyll[d] <- sd(log10(my.df$Avg_Chloro), na.rm = T)
  phyto.df$volatility.fcm.AF[d] <- sd(my.df$total_ml_vol_corrected_AF, na.rm = T)
  # phyto.df$volatility.O2bio.estimated[d] <- sd(my.df$O2bio.estimated, na.rm = T)
  phyto.df$volatility.sum.count[d] <- sd(my.df$sum.phyto.count, na.rm = T)
  phyto.df$volatility.log10.sum.count[d] <- sd(log10(my.df$sum.phyto.counts), na.rm = T)
  phyto.df$volatility.temperature[d] <- sd(my.df$SURF_TEMP_C, na.rm = T)
  phyto.df$volatility.temperature.anomaly[d] <- sd(my.df$SURF_TEMP_ANOM, na.rm = T)
  phyto.df$volatility.richness[d] <- sd(my.df$richness, na.rm = T)
  phyto.df$volatility.shannon.div[d] <- sd(my.df$shannon.div, na.rm = T)
  phyto.df$volatility.nitrate[d] <- sd(my.df$Nitrate, na.rm = T)
  phyto.df$volatility.phosphate[d] <- sd(my.df$Phosphate, na.rm = T)
  phyto.df$volatility.silicate[d] <- sd(my.df$Silicate, na.rm = T)
  phyto.df$volatility.beuti[d] <- sd(my.df$BEUTI.33N, na.rm = T)
  
  
}

phyto.df.o2biodates$volatility.O2bio.estimated <- NA

for(d in 1:nrow(phyto.df.o2biodates)){
  
  my.date <- phyto.df.o2biodates$Date[d]
  
  my.df <- phyto.df.o2biodates[which(phyto.df.o2biodates$Date >= (my.date - 60*60*24*15) & phyto.df.o2biodates$Date <= (my.date + 60*60*24*15)),]
  
  phyto.df.o2biodates$volatility.O2bio.estimated[d] <- sd(my.df$O2bio.estimated, na.rm = T)
  
}


# ---- plot volatility ----


# plot.vol.O2bio.predalt <- ggplot(data = phyto.df) +
#   geom_line(aes(x = Date, y = volatility.O2bio.predicted.altered)) +
#   labs(y = "Predicted [O2]bio Volatility (uM)") +
#   theme_bw()

# plot.vol.O2bio.est <- ggplot(data = phyto.df.o2biodates) +
#   geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
#   geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
#   geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
#   geom_line(aes(x = Date, y = volatility.O2bio.estimated), color = "black", lwd = 1) +
#   labs(y = "[O2]bio \n Volatility (uM)") +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))
# 
# plot.vol.O2bio.est.box <- ggplot(phyto.df.o2biodates) +
#   geom_boxplot(aes(x = Blob.status, y = volatility.O2bio.estimated), fill = "black") +
#   scale_fill_manual(values = my.colors.blob.2022) +
#   labs(x = "Time Period", y = "[O2]bio \n Volatility (uM)") +
#   guides(fill = "none") +
#   # geom_text(
#   #   data = label_df,
#   #   aes(x = Blob.status, y = (y_pos), label = Letter),
#   #   inherit.aes = FALSE) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.O2bio.predalt <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = volatility.O2bio.predicted.altered), color = "black", lwd = 1) +
  labs(y = "Predicted \n [O2]bio \n Volatility (uM)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.O2bio.predalt.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.O2bio.predicted.altered), fill = "black") +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Predicted \n [O2]bio \n Volatility (uM)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


plot.vol.beuti <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = volatility.beuti), color = "grey30", lwd = 1) +
  labs(y = "BEUTI Index \n 33N Volatility \n (mmol/m/sec)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.beuti.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.beuti), fill = "grey30") +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "BEUTI Index \n 33N Volatility \n (mmol/m/sec)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


# plot.vol.chlorophyll <- ggplot(data = phyto.df) +
#   geom_line(aes(x = Date, y = volatility.chlorophyll)) +
#   theme_bw()

plot.vol.log10chlorophyll <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = volatility.log10chlorophyll), color = my.colorblind.colors[3], lwd = 1) +
  labs(y = "log10 \n Chlorophyll \n Volatility (ug/L)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.chla.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.log10chlorophyll), fill = my.colorblind.colors[3]) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "log10 \n Chlorophyll \n Volatility (ug/L)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

# plot.vol.FCM <- ggplot(data = phyto.df) +
#   geom_line(aes(x = Date, y = volatility.fcm.AF)) +
#   theme_bw()

# plot.vol.sum.counts <- ggplot(data = phyto.df) +
#   geom_line(aes(x = Date, y = volatility.sum.count)) +
#   theme_bw()

plot.vol.log10.sum.counts <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = volatility.log10.sum.count), color = my.colorblind.colors[5], lwd = 1) +
  labs(y = "log10 Total \n Phytoplankton Count \n Volatility") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.log10.sum.counts.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.log10.sum.count), fill = my.colorblind.colors[5]) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "log10 Total \n Phytoplankton Count \n Volatility") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.temp <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = volatility.temperature), color = my.colorblind.colors[2], lwd = 1) +
  labs(y = "Surface Water \n Temperature \n Volatility (C)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.temp.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.temperature), fill = my.colorblind.colors[2]) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Surface Water \n Temperature \n Volatility (C)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.temp.anomaly <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = volatility.temperature.anomaly), color = my.colorblind.colors[4], lwd = 1) +
  labs(y = "Surface Water \n Temperature Anomaly \n Volatility (C)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.temp.anomaly.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.temperature.anomaly), fill = my.colorblind.colors[4]) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Surface Water \n Temperature Anomaly \n Volatility (C)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.richness <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = volatility.richness), color = my.colorblind.colors[6], lwd = 1) +
  labs(y = "Phytoplankton \n Richness \n Volatility") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.richness.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.richness), fill = my.colorblind.colors[6]) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Phytoplankton \n Richness \n Volatility") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.shannondiv <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = volatility.shannon.div), color = my.colorblind.colors[1], lwd = 1) +
  labs(y = "Phytoplankton \n Diversity \n Volatility") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.shannondiv.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.shannon.div), fill = my.colorblind.colors[1]) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Phytoplankton \n Diversity \n Volatility") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.nitrate <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = volatility.nitrate), color = "grey69", lwd = 1) +
  labs(y = "Nitrate \n Volatility (uM)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.nitrate.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.nitrate), fill = "white") +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Nitrate \n Volatility (uM)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.phosphate <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = volatility.phosphate), color = "grey69", lwd = 1) +
  labs(y = "Phosphate \n Volatility (uM)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.phosphate.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.phosphate), fill = "white") +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Phosphate \n Volatility (uM)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.silicate <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = volatility.silicate), color = "grey69", lwd = 1) +
  labs(y = "Silicate \n Volatility (uM)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.vol.silicate.box <- ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.silicate), fill = "white") +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Silicate \n Volatility (uM)") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = (y_pos), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

# plot.vol.O2bio.est + plot.vol.log10.sum.counts + plot.vol.log10chlorophyll +plot.vol.temp.anomaly + plot.vol.temp + plot.vol.richness + plot.vol.shannondiv + plot_layout(ncol = 1)



plot.temp <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = SURF_TEMP_C), color = my.colorblind.colors[2], lwd = 1) +
  labs(y = "Surface Water \n Temperature (C)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.temp.anomaly <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = SURF_TEMP_ANOM), color = my.colorblind.colors[4], lwd = 1) +
  labs(y = "Surface Water \n Temperature \n Anomaly (C)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.log10chlorophyll <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = log10(Avg_Chloro)), color = my.colorblind.colors[3], lwd = 1) +
  labs(y = "log10 \n Chlorophyll \n (ug/L)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.log10.sum.counts <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = log10(sum.phyto.counts)), color = my.colorblind.colors[5], lwd = 1) +
  labs(y = "log10 Total \n Phytoplankton \n Count") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

# plot.O2bio.est <- ggplot(data = phyto.df.o2biodates) +
#   geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
#   geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
#   geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
#   geom_line(aes(x = Date, y = O2bio.estimated), color = "black", lwd = 1) +
#   labs(y = "[O2]bio (uM)") +
#   scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
#   theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.O2bio.predalt <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = O2bio.predicted.altered), color = "black", lwd = 1) +
  labs(y = "Predicted \n [O2]bio (uM)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.beuti <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = BEUTI.33N), color = "grey30", lwd = 1) +
  labs(y = "BEUTI Index \n 33N \n (mmol/m/sec)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))


plot.richness <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = richness), color = my.colorblind.colors[6], lwd = 1) +
  labs(y = "Phytoplankton \n Richness") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.shannondiv <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = shannon.div), color = my.colorblind.colors[1], lwd = 1) +
  labs(y = "Phytoplankton \n Diversity") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.nitrate <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = Nitrate), color = "grey69", lwd = 1) +
  labs(y = "Nitrate (uM)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.phosphate <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = Phosphate), color = "grey69", lwd = 1) +
  labs(y = "Phosphate (uM)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))  +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

plot.silicate <- ggplot(data = phyto.df) +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.01) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = Date, y = Silicate), color = "grey69", lwd = 1) +
  labs(y = "Silicate (uM)") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))

# plot.O2bio.est + plot.log10.sum.counts + plot.log10chlorophyll + plot.temp.anomaly + plot.temp + plot.richness + plot.shannondiv + plot_layout(ncol = 1)


# ----- plot big ole plot ----

## BIG PLOT!! 
big.ole.plot <- plot.o2bio.predalt.box + plot.O2bio.predalt + plot.vol.O2bio.predalt + plot.vol.O2bio.predalt.box + 
  plot.richness.box + plot.richness +  plot.vol.richness + plot.vol.richness.box +
  plot.shannon.div.box + plot.shannondiv + plot.vol.shannondiv + plot.vol.shannondiv.box +
  plot.phytocount.box + plot.log10.sum.counts + plot.vol.log10.sum.counts + plot.vol.log10.sum.counts.box + 
  plot.chla.box + plot.log10chlorophyll + plot.vol.log10chlorophyll + plot.vol.chla.box + 
  plot.tempanom.box + plot.temp.anomaly + plot.vol.temp.anomaly + plot.vol.temp.anomaly.box + 
  plot.temp.box + plot.temp + plot.vol.temp + plot.vol.temp.box + 
  plot.nitrate.box + plot.nitrate + plot.vol.nitrate + plot.vol.nitrate.box +
  plot.phosphate.box + plot.phosphate + plot.vol.phosphate + plot.vol.phosphate.box +
  plot.silicate.box + plot.silicate + plot.vol.silicate + plot.vol.silicate.box +
  plot.beuti.box + plot.beuti + plot.vol.beuti + plot.vol.beuti.box +
  plot_layout(ncol = 4, axes = "collect", widths = c(1,3,3,1))

ggsave("Figures/2026-04-19_big_ole_plot.pdf", width = 20, height = 20, units = "in")





summary(lm(phyto.df$volatility.log10chlorophyll~phyto.df$volatility.log10.sum.count))
ggplot(data = phyto.df) +
  geom_point(aes(x = volatility.log10.sum.count, y = volatility.log10chlorophyll)) +
  geom_smooth(aes(x = volatility.log10.sum.count, y = volatility.log10chlorophyll), method = "lm", se = F) +
  theme_bw()



my.colors.blob.2022 <- my.colorblind.colors[c(3,2,4,6)]


## Chlorophyll volatility ##
my.aov <- aov(phyto.df$volatility.log10chlorophyll~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$volatility.log10chlorophyll, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.log10chlorophyll, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Volatility log10 Chlorophyll") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = log10(y_pos)*2, label = Letter),
    inherit.aes = FALSE) +
  theme_bw()



## Phyto counts volatility ##
my.aov <- aov(phyto.df$volatility.log10.sum.count~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$volatility.log10.sum.count, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.log10.sum.count, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Volatility log10 Counts") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = log10(y_pos)*2, label = Letter),
    inherit.aes = FALSE) +
  theme_bw()



## Richness volatility ##
my.aov <- aov(phyto.df$volatility.richness~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$volatility.richness, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.richness, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Volatility Richness") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = y_pos*2, label = Letter),
    inherit.aes = FALSE) +
  theme_bw()



## Shannon diversity volatility ##
my.aov <- aov(phyto.df$volatility.shannon.div~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$volatility.shannon.div, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.shannon.div, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Volatility Diversity") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = y_pos*2, label = Letter),
    inherit.aes = FALSE) +
  theme_bw()


## water temp volatility ##
my.aov <- aov(phyto.df$volatility.temperature~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$volatility.temperature, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.temperature, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Volatility Water Temperature") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = y_pos*2, label = Letter),
    inherit.aes = FALSE) +
  theme_bw()


## water temp anoamly volatility ##
my.aov <- aov(phyto.df$volatility.temperature.anomaly~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$volatility.temperature.anomaly, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.temperature.anomaly, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Volatility Water Temperature Anomaly") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = y_pos*2, label = Letter),
    inherit.aes = FALSE) +
  theme_bw()

## nitrate volatility ##
my.aov <- aov(phyto.df$volatility.nitrate~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$volatility.nitrate, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.nitrate, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Volatility Nitrate") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = y_pos*2, label = Letter),
    inherit.aes = FALSE) +
  theme_bw()


## phosphate volatility ##
my.aov <- aov(phyto.df$volatility.phosphate~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$volatility.phosphate, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.phosphate, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Volatility Phosphate") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = y_pos*2, label = Letter),
    inherit.aes = FALSE) +
  theme_bw()


## silicate volatility ##
my.aov <- aov(phyto.df$volatility.silicate~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$volatility.silicate, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.silicate, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Volatility Silicate") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = y_pos*2, label = Letter),
    inherit.aes = FALSE) +
  theme_bw()


## predicted O2bio volatility ##
my.aov <- aov(phyto.df$volatility.O2bio.predicted.altered~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$volatility.O2bio.predicted.altered, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.O2bio.predicted.altered, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Volatility Predicted [O2]bio") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = y_pos*2, label = Letter),
    inherit.aes = FALSE) +
  theme_bw()




# ---- COMPARISONS ----
# summary(lm(phyto.df$volatility.O2bio.predicted.altered~phyto.df$volatility.log10chlorophyll))
# ggplot(phyto.df) +
#   geom_point(aes(x = volatility.log10chlorophyll, y = volatility.O2bio.predicted.altered)) +
#   geom_smooth(aes(x = volatility.log10chlorophyll, y = volatility.O2bio.predicted.altered), method = "lm", se = F) +
#   labs(x = "Volatility log10 Chlorophyll", y = "Volatility Predicted [O2]bio") +
#   guides(fill = "none") +
#   # geom_text(
#   #   data = label_df,
#   #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
#   #   inherit.aes = FALSE) +
#   theme_bw()

# summary(lm(phyto.df$volatility.O2bio.estimated~phyto.df$volatility.log10chlorophyll))
# ggplot(phyto.df) +
#   geom_point(aes(x = volatility.log10chlorophyll, y = volatility.O2bio.estimated)) +
#   geom_smooth(aes(x = volatility.log10chlorophyll, y = volatility.O2bio.estimated), method = "lm", se = F) +
#   labs(x = "Volatility log10 Chlorophyll", y = "Volatility Estimated [O2]bio - Model 1") +
#   guides(fill = "none") +
#   # geom_text(
#   #   data = label_df,
#   #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
#   #   inherit.aes = FALSE) +
#   theme_bw()
# 
# summary(lm(phyto.df$volatility.O2bio.estimated~phyto.df$volatility.fcm.AF))
# ggplot(phyto.df) +
#   geom_point(aes(x = volatility.fcm.AF, y = volatility.O2bio.estimated)) +
#   geom_smooth(aes(x = volatility.fcm.AF, y = volatility.O2bio.estimated), method = "lm", se = F) +
#   labs(x = "Volatility log10 Chlorophyll", y = "Volatility Estimated [O2]bio - Model 1") +
#   guides(fill = "none") +
#   # geom_text(
#   #   data = label_df,
#   #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
#   #   inherit.aes = FALSE) +
#   theme_bw()


phyto.df <- merge(phyto.df, phyto.nmds, by = "Date", all.x = T, all.y = F)

ggplot() +
  geom_point(data = phyto.df[which(is.na(phyto.df$O2bio.estimated) == T),], aes(x = MDS1, y = MDS2), color = "grey69", alpha = 0.3, size = 3) +
  geom_point(data = phyto.df[which(is.na(phyto.df$O2bio.estimated) == F),], aes(x = MDS1, y = MDS2, color = O2bio.estimated), alpha = 1, size = 3) +
  # scale_color_manual(values = rainbow(12)) +
  # scale_fill_manual(values = rainbow(12)) +
  scale_color_viridis_c() +
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "O2bio.estimated",  x = "Dim 1", y = "Dim 2") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 

ggplot() +
  geom_point(data = phyto.df[which(is.na(phyto.df$volatility.O2bio.estimated) == T),], aes(x = MDS2, y = MDS3), color = "grey69", alpha = 0.3, size = 3) +
  geom_point(data = phyto.df[which(is.na(phyto.df$volatility.O2bio.estimated) == F),], aes(x = MDS2, y = MDS3, color = volatility.O2bio.estimated), alpha = 1, size = 3) +
  # scale_color_manual(values = rainbow(12)) +
  # scale_fill_manual(values = rainbow(12)) +
  scale_color_viridis_c() +
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Volatility \n O2bio.estimated",  x = "Dim 2", y = "Dim 3") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 

summary(lm(phyto.df$volatility.O2bio.estimated~phyto.df$MDS1))
summary(lm(phyto.df$volatility.O2bio.estimated~phyto.df$MDS2))
summary(lm(phyto.df$volatility.O2bio.estimated~phyto.df$MDS3))
summary(lm(phyto.df$volatility.O2bio.estimated~phyto.df$MDS4))

ggplot(phyto.df) +
  geom_point(aes(x = MDS2, y = volatility.O2bio.estimated)) +
  geom_smooth(aes(x = MDS2, y = volatility.O2bio.estimated), method = "lm", se = F) +
  labs(x = "Dim 2", y = "Volatility Estimated [O2]bio - Model 1") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw()

ggplot(phyto.df) +
  geom_point(aes(x = MDS3, y = volatility.O2bio.estimated)) +
  geom_smooth(aes(x = MDS3, y = volatility.O2bio.estimated), method = "lm", se = F) +
  labs(x = "Dim 3", y = "Volatility Estimated [O2]bio - Model 1") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw()


summary(lm(phyto.df$O2bio.estimated~phyto.df$MDS1))
summary(lm(phyto.df$O2bio.estimated~phyto.df$MDS2))
summary(lm(phyto.df$O2bio.estimated~phyto.df$MDS3))
summary(lm(phyto.df$O2bio.estimated~phyto.df$MDS4))

ggplot() +
  geom_point(data = phyto.df, aes(x = MDS2, y = MDS3, color = factor(Year)), alpha = 0.3, size = 3) +
  # scale_color_manual(values = rainbow(12)) +
  # scale_fill_manual(values = rainbow(12)) +
  # scale_color_viridis_c() +
  scale_color_viridis_d() +
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Year",  x = "Dim 2", y = "Dim 3") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 



## sum phyto counts vs chlorophyll

summary(lm(phyto.df$Avg_Chloro~phyto.df$sum.phyto.counts))
ggplot(phyto.df) +
  geom_point(aes(x = Avg_Chloro, y = sum.phyto.counts)) +
  geom_smooth(aes(x = Avg_Chloro, y = sum.phyto.counts), method = "lm", se = F) +
  labs(x = "Dim 2", y = "Volatility Estimated [O2]bio - Model 1") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw()

summary(lm(log10(phyto.df$Avg_Chloro)~log10(phyto.df$sum.phyto.counts)))
ggplot(phyto.df) +
  geom_point(aes(x = log10(Avg_Chloro), y = log10(sum.phyto.counts))) +
  geom_smooth(aes(x = log10(Avg_Chloro), y = log10(sum.phyto.counts)), method = "lm", se = F) +
  labs(x = "Dim 2", y = "Volatility Estimated [O2]bio - Model 1") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw()



# ---- investigate mechanism for volatility ----

# summary(lm(phyto.df$volatility.log10chlorophyll~phyto.df$SURF_TEMP_C))
# ggplot(data = phyto.df) +
#   geom_point(aes(x = SURF_TEMP_C, y = volatility.log10chlorophyll)) +
#   geom_smooth(aes(x = SURF_TEMP_C, y = volatility.log10chlorophyll), method = "lm", se = F) +
#   theme_bw()
# 
# summary(lm(phyto.df$volatility.log10chlorophyll~phyto.df$volatility.temperature))
# ggplot(data = phyto.df) +
#   geom_point(aes(x = volatility.temperature, y = volatility.log10chlorophyll)) +
#   geom_smooth(aes(x = volatility.temperature, y = volatility.log10chlorophyll), method = "lm", se = F) +
#   theme_bw()


my.cor.variables <- c("O2bio.predicted.altered", "richness", "shannon.div", "sum.phyto.counts", "Avg_Chloro", "SURF_TEMP_C", "SURF_TEMP_ANOM", "Nitrate", "Phosphate", "Silicate", "BEUTI.33N", colnames(phyto.df)[which(substr(colnames(phyto.df),1,10) == "volatility")])
my.cor.variables <- my.cor.variables[which(my.cor.variables != "volatility.fcm.AF")]
my.cor.variables <- my.cor.variables[which(my.cor.variables != "volatility.chlorophyll")]
my.cor.variables <- my.cor.variables[which(my.cor.variables != "volatility.sum.count")]


my.cor.df <- phyto.df[,which(colnames(phyto.df) %in% my.cor.variables)]

my.cor.df$log10.sum.count <- log10(my.cor.df$sum.phyto.counts)
my.cor.df$log10.chlorophyll <- log10(my.cor.df$Avg_Chloro)
my.cor.df <- my.cor.df[,-which(colnames(my.cor.df) %in% c("sum.phyto.counts", "Avg_Chloro"))]


my.cor.df.reordered <- my.cor.df[, c("O2bio.predicted.altered", "richness", "shannon.div", "log10.sum.count", "log10.chlorophyll", "SURF_TEMP_ANOM", "SURF_TEMP_C", "Nitrate", "Phosphate", "Silicate", "BEUTI.33N", "volatility.O2bio.predicted.altered", "volatility.richness", "volatility.shannon.div", "volatility.log10.sum.count", "volatility.log10chlorophyll", "volatility.temperature.anomaly", "volatility.temperature", "volatility.nitrate", "volatility.phosphate", "volatility.silicate", "volatility.beuti")]

colnames(my.cor.df.reordered) <- c("Predicted [O2]bio", "Richness", "Shannon Diversity", "log10 Phyto Count", "log10 Chlorophyll", "Surface Water Temperature Anomaly", "Surface Water Temperature", "Nitrate", "Phosphate", "Silicate", "BEUTI Index 33N", "Predicted [O2]bio Volatility", "Richness Volatility", "Shannon Diversity Volatility", "log10 Phyto Count Volatility", "log10 Chlorophyll Volatility",  "Surface Water Temperature Anomaly Volatility", "Surface Water Temperature Volatility", "Nitrate Volatility", "Phosphate Volatility", "Silicate Volatility", "BEUTI Index 33N Volatility")

my.cor.results <- cor(my.cor.df.reordered, method = "pearson", use = "pairwise.complete.obs")




my.tl.colors <- c("black", my.colorblind.colors[6], my.colorblind.colors[1], my.colorblind.colors[5], my.colorblind.colors[3], my.colorblind.colors[4], my.colorblind.colors[2], "grey69", "grey69", "grey69", "grey30")

par(font = 2)
corrplot(my.cor.results, diag = F, method = "square", tl.cex = 0.5, mar = c(0,0,0,0), col = viridis(100), rect.col = "black", tl.col = my.tl.colors, cl.)

my.cor.results <- as.data.frame(my.cor.results)
my.cor.results$X <- rownames(my.cor.results)
my.cor.results.long <- my.cor.results %>% pivot_longer(cols = colnames(my.cor.results)[-which(colnames(my.cor.results) == "X")], names_to = "Y", values_to = "Perason.r")




my.cor.results.list <- rcorr(as.matrix(my.cor.df.reordered), type = "pearson")
my.cor.results.p <- my.cor.results.list$P
my.cor.results.p <- as.data.frame(my.cor.results.p)
my.cor.results.p$X <- rownames(my.cor.results.p)
my.cor.results.p.long <- my.cor.results.p %>% pivot_longer(cols = colnames(my.cor.results.p)[-which(colnames(my.cor.results.p) == "X")], names_to = "Y", values_to = "Perason.p")


my.cor.results.full <- merge(my.cor.results.long, my.cor.results.p.long, by = c("X","Y"))
my.cor.results.full <- my.cor.results.full[which(my.cor.results.full$Perason.p < 0.05),]
my.cor.results.full$abs.Pearson.r <- abs(my.cor.results.full$Perason.r)
my.cor.results.full$R2 <- (my.cor.results.full$Perason.r)^2



# ---- jaccard community turnover -----

my.presence.absence.mat <- (phyto.mat > 0) + 0

# Compute turnover (replacement)
pair_comp <- beta.pair(my.presence.absence.mat, index.family="jaccard")

# Extract the turnover component
turnover_jaccard <- pair_comp$beta.jtu
print(turnover_jaccard)

try.it <- as.matrix(turnover_jaccard)

super_diag <- diag(try.it[-nrow(try.it), -1])

ggplot() +
  geom_rect(aes(xmin=parse_date_time("2020-03-20", orders = "Ymd"), xmax = parse_date_time("2020-06-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'springgreen', alpha = 0.5) +
  geom_rect(aes(xmin=parse_date_time("2013-10-01", orders = "Ymd"), xmax = parse_date_time("2016-04-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_rect(aes(xmin=parse_date_time("2014-01-01", orders = "Ymd"), xmax = parse_date_time("2016-01-01", orders = "ymd"), ymin = -Inf, ymax = Inf), fill = 'pink', alpha = 0.5) +
  geom_line(aes(x = parse_date_time(rownames(try.it), orders = "Ymd"), y = c(NA,super_diag)), lwd = 1) +
  labs(y = "Jaccard Turnover Index", x = "Date") +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 


# ---- taxa variability between time periods ----

my.abund.taxa <- colnames(phyto.mat)


phytos.combo$Blob.status[which(phytos.combo$Year >= 2022)] <- "Post2022Shift"
phytos.combo$Blob.status <- factor(phytos.combo$Blob.status, levels = c("PreBlob", "Blob", "PostBlob", "Post2022Shift"))


for(t in 1:length(my.abund.taxa)){
  
  my.taxon <- my.abund.taxa[t]
  
  my.abund <- phytos.combo[,which(colnames(phytos.combo) == my.taxon)]
  # my.log10.abund[which(my.log10.abund < 0)] <- 0
  
  my.aov <- aov(my.abund~phytos.combo$Blob.status)
  
  my.tukey <- TukeyHSD(my.aov)
  
  my.tukey.results <- my.tukey$`phytos.combo$Blob.status`
  
  my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
  my.df$period.comparison <- rownames(my.df)
  
  my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters
  
  y_max <- tapply(my.abund, phytos.combo$Blob.status, max, na.rm = TRUE)
  
  label_df <- data.frame(
    Blob.status = names(my.letters),
    Letter = my.letters,
    y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
  )
  
  if(nrow(my.df) > 0){
    
    png(filename = paste0("Figures/taxa_shift_boxplots/", "2026-04-09_taxa_shift_boxplot_", my.taxon, ".png"))
    
    c <- ggplot(data = phytos.combo, 
                aes(x = Blob.status, y = log10(phytos.combo[,which(colnames(phytos.combo) == my.taxon)]), fill = Blob.status, group = Blob.status)) +
      geom_boxplot() +
      # geom_signif(comparisons = list(c("PreBlob", "Blob"),
      #                                c("Blob", "PostBlob"),
      #                                c("PostBlob", "Post2022Shift")),
      #             map_signif_level = T) +
      geom_text(
        data = label_df,
        aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
        inherit.aes = FALSE) +
      scale_fill_manual(values = my.colors.blob.2022) +
      theme_bw() +
      labs(x = "Time Period", y = "log10(Phyto Count)") +
      theme(legend.position = "none") +
      ggtitle(my.taxon)
    
    print(c)
    
    dev.off()
    
  }
  
  
  
}


my.taxon <- "O2bio.predicted.altered"

my.O2bio <- phytos.combo[,which(colnames(phytos.combo) == my.taxon)]

my.aov <- aov(my.O2bio~phytos.combo$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phytos.combo$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(my.O2bio, phytos.combo$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(data = phytos.combo, 
       aes(x = Blob.status, y = my.O2bio, fill = Blob.status, group = Blob.status)) +
  geom_boxplot() +
  # geom_signif(comparisons = list(c("PreBlob", "Blob"),
  #                                c("Blob", "PostBlob"),
  #                                c("PostBlob", "Post2022Shift")),
  #             map_signif_level = T) +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = y_pos, label = Letter),
    inherit.aes = FALSE) +
  scale_fill_manual(values = my.colors.blob.2022) +
  theme_bw() +
  labs(x = "Time Period", y = "Predicted O2bio") +
  theme(legend.position = "none") +
  ggtitle(my.taxon)


saveRDS(phytos.combo, file = "2026-04-09_phytos_combo_nmds_o2bio.rds")


# ---- differential abundance of taxa between time periods ----

## lots of this code is from 2023-11-06_deseq(1).R


combo.full.mat <- phytos.combo[,which(colnames(phytos.combo) %in% predictors)]

rownames(combo.full.mat) <- phytos.combo$Date

## Make test data if you forget which way is which
# combo.full.mat$TEST <- 0
# combo.full.mat$TEST[which(as.numeric(substr(rownames(combo.full.mat), start = 1, stop = 4)) < 2022)] <- 10000

#We will convert our table to DESeqDataSet object
countData = round(as(combo.full.mat, "matrix"), digits = 0)
countData[is.na(countData)] <- 0
# countData <- countData[,which(colSums(countData) >= 100)]

## We will add 1 to the countData otherwise DESeq will fail with the error:
## estimating size factors
## Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  :
## every gene contains at least one zero, cannot compute log geometric means
countData<-(t(countData+1)) 


# metadata <- as.data.frame(combo.full$Dataset)
# colnames(metadata) <- "Dataset"


metadata <- data.frame(Date = parse_date_time(rownames(combo.full.mat), orders = "Ymd"))
metadata$Blob.status <- "Blob"
metadata$Blob.status[which(year(metadata$Date) <= 2013)] <- "PreBlob"
metadata$Blob.status[which(year(metadata$Date) >= 2016)] <- "PostBlob"
metadata$Blob.status[which(year(metadata$Date) >= 2022)] <- "Post2022Shift"
my.rownames <- metadata$Date
metadata <- as.data.frame(metadata[,-1])
rownames(metadata) <- my.rownames
colnames(metadata) <- "Blob.status"

# metadata$Dataset[which(substr(metadata$Dataset, 1,1) == "W")] <- "WEAllen"
# metadata$heatwave.perc.90 <- factor(metadata$Dataset)
# rownames(metadata) <- combo.full$Date


## subset to shift of interest here ##
metadata$Blob.status <- factor(metadata$Blob.status, levels = c("PreBlob", "Blob", "PostBlob", "Post2022Shift"))
# index <- which(metadata$Blob.status %in% c("PreBlob", "Blob"))
# index <- which(metadata$Blob.status %in% c("Blob", "PostBlob"))
index <- which(metadata$Blob.status %in% c("PostBlob", "Post2022Shift"))

countData <- countData[,index]
meta.dates <- rownames(metadata)
metadata <- as.data.frame(metadata[index,])
rownames(metadata) <- meta.dates[index]
colnames(metadata) <- "Blob.status"


dds <- DESeqDataSetFromMatrix(countData, metadata, as.formula(~Blob.status))
## ignore warning about characters/factors

#Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
data_deseq_test = DESeq(dds)

## Extract the results
res = results(data_deseq_test, cooksCutoff = FALSE)
res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]), OTU = rownames(res))

sig =0.05  
fold = 0
plot.point.size = 2
label=T
tax.display = NULL
tax.aggregate = "OTU"

res_tax_sig = subset(res_tax, padj < sig & fold < abs(log2FoldChange))

res_tax_sig <- res_tax_sig[order(res_tax_sig$padj),]

## Plot the data MA plot

### MA plot
res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
res_tax$Significant[is.na(res_tax$Significant)] <- "No"
p1 <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) +
  geom_point(size = plot.point.size) +
  scale_x_log10() +
  scale_color_manual(values=c("black", "red")) +
  labs(x = "Mean abundance", y = "Log2 fold change")+theme_bw()

p1

#Running theif statement plots for names that are significant 
if(label == T){
  if (!is.null(tax.display)){
    rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
  } else {
    rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
  }
  p1 <- p1 + geom_text(data = subset(rlab, Significant == "Yes"), aes(label = Display), size = 4, vjust = 1)
}

p1


res_tax_sig <- res_tax_sig[order(res_tax_sig$log2FoldChange, decreasing = F),]
res_tax_sig$Predictors <- factor(res_tax_sig$OTU, levels = res_tax_sig$OTU[order(res_tax_sig$log2FoldChange, decreasing = T)])

res_tax_sig$Group <- NA
res_tax_sig$Group[which(res_tax_sig$Predictors %in% unique.diatoms)] <- "Diatom"
res_tax_sig$Group[which(res_tax_sig$Predictors %in% unique.dinos)] <- "Dinoflagellate"


ggplot(data = res_tax_sig) +
  # annotate(geom = "rect", xmin = min(as.numeric(res_tax_sig$Predictors))-0.5, xmax = max(as.numeric(res_tax_sig$Predictors))+0.5, ymin = -Inf, ymax = 0, fill = my.colorblind.colors[3], alpha = 0.4) +
  # annotate(geom = "rect", xmin = min(as.numeric(res_tax_sig$Predictors))-0.5, xmax = max(as.numeric(res_tax_sig$Predictors))+0.5, ymin = 0, ymax = Inf, fill = my.colorblind.colors[2], alpha = 0.4) +
  # annotate(geom = "rect", xmin = min(as.numeric(res_tax_sig$Predictors))-0.5, xmax = max(as.numeric(res_tax_sig$Predictors))+0.5, ymin = -Inf, ymax = 0, fill = my.colorblind.colors[2], alpha = 0.4) +
  # annotate(geom = "rect", xmin = min(as.numeric(res_tax_sig$Predictors))-0.5, xmax = max(as.numeric(res_tax_sig$Predictors))+0.5, ymin = 0, ymax = Inf, fill = my.colorblind.colors[4], alpha = 0.4) +
  annotate(geom = "rect", xmin = min(as.numeric(res_tax_sig$Predictors))-0.5, xmax = max(as.numeric(res_tax_sig$Predictors))+0.5, ymin = -Inf, ymax = 0, fill = my.colorblind.colors[4], alpha = 0.4) +
  annotate(geom = "rect", xmin = min(as.numeric(res_tax_sig$Predictors))-0.5, xmax = max(as.numeric(res_tax_sig$Predictors))+0.5, ymin = 0, ymax = Inf, fill = my.colorblind.colors[6], alpha = 0.4) +
  geom_bar(aes(x = Predictors, y = log2FoldChange, fill = Group), stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("deepskyblue2", "deepskyblue4")) +
  coord_flip() +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  # ggtitle(label = "Blob Shift") + theme(plot.title = element_text(size = 14, face = "bold", hjust = 1)) +
  # labs(tag = "A", x = "Taxon", y = "log2 Fold Change") + theme(plot.tag = element_text(size = 18, face = "bold"))
  # ggtitle(label = "Post-Blob Shift") + theme(plot.title = element_text(size = 14, face = "bold", hjust = 1)) +
  # labs(tag = "B", x = "Taxon", y = "log2 Fold Change") + theme(plot.tag = element_text(size = 18, face = "bold"))
  ggtitle(label = "2022 Shift") + theme(plot.title = element_text(size = 14, face = "bold", hjust = 1)) +
  labs(tag = "C", x = "Taxon", y = "log2 Fold Change") + theme(plot.tag = element_text(size = 18, face = "bold"))



length(which(res_tax_sig$log2FoldChange > 0 & res_tax_sig$Group == "Diatom"))
length(which(res_tax_sig$log2FoldChange > 0))
length(which(res_tax_sig$log2FoldChange > 0 & res_tax_sig$Group == "Diatom"))/length(which(res_tax_sig$log2FoldChange > 0))

length(which(res_tax_sig$log2FoldChange < 0 & res_tax_sig$Group == "Diatom"))
length(which(res_tax_sig$log2FoldChange < 0))
length(which(res_tax_sig$log2FoldChange < 0 & res_tax_sig$Group == "Diatom"))/length(which(res_tax_sig$log2FoldChange > 0))


as.numeric(res_tax_sig$Predictors)










