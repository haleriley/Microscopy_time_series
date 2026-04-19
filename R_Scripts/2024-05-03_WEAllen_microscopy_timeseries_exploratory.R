# exploring W. E. Allen phytoplankton microscopy time series
# RJH
# 2024-05-03

# ---- library ----

library(tidyverse)
library(lubridate)
library(janitor)
library(vegan)


# ---- read in data ----

setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/")

diatoms <- read.csv("AllenData/SIOdiat.csv")
dinos <- read.csv("AllenData/SIOdino.csv", check.names = F)

# ---- clean and reformat data ----

diatoms <- clean_names(diatoms)
dinos <- clean_names(dinos)

diatoms$date <- parse_date_time(paste(diatoms$year, diatoms$month, diatoms$day), orders = "Ymd")
dinos$date <- parse_date_time(paste(dinos$year, dinos$month, dinos$day), orders = "Ymd")

diatoms[is.na(diatoms)] <- 0
dinos[is.na(dinos)] <- 0

dinos <- dinos[,c(ncol(dinos), 1:ncol(dinos)-1)]
diatoms <- diatoms[,c(ncol(diatoms), 1:ncol(diatoms)-1)]

dinos <- dinos[which(rowSums(dinos[5:ncol(dinos)]) >= 10), which(colSums(dinos[5:ncol(dinos)]) >= 10)]
dinos.long <- dinos %>% pivot_longer(cols = colnames(dinos)[5:ncol(dinos)], names_to = "taxon", values_to = "count")

diatoms <- diatoms[which(rowSums(diatoms[5:ncol(diatoms)]) >= 10), which(colSums(diatoms[5:ncol(diatoms)]) >= 10)]
diatoms.long <- diatoms %>% pivot_longer(cols = colnames(diatoms)[5:ncol(diatoms)], names_to = "taxon", values_to = "count")


dinos.long <- dinos.long %>% separate(taxon, c("genus", "species")) %>% group_by(genus, date) %>% summarize(count = mean(count))

ggplot(data = dinos.long) +
  geom_line(aes(x = date, y = count, color = genus), size = 1) +
  theme_bw() +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Date", y = "Count", color = "Genus") +
  ggtitle("W. E. Allen Phytoplankton Microscopy Time Series - Dinoflagellates") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  theme(panel.grid.minor.x = element_blank())

diatoms.long <- diatoms.long %>% separate(taxon, c("genus", "species")) %>% group_by(genus, date) %>% summarize(count = mean(count))

ggplot(data = diatoms.long) +
  geom_line(aes(x = date, y = count, color = genus), size = 1) +
  theme_bw() +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Date", y = "Count", color = "Genus") +
  ggtitle("W. E. Allen Phytoplankton Microscopy Time Series - Diatoms") +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) +
  theme(panel.grid.minor.x = element_blank())



phytos <- merge(diatoms, dinos, by = c("date", "year", "month", "day"), all = F)

phytos.long <- phytos %>% pivot_longer(cols = colnames(phytos)[5:ncol(phytos)], names_to = "taxon", values_to = "count")

phytos <- diatoms

my.rownames <- phytos$date
phytos <- phytos[,-c(1:3, ncol(phytos))]
rownames(phytos) <- my.rownames

phytos <- phytos[which(rowSums(phytos) != 0),]
phytos <- phytos[,which(colSums(phytos) != 0)]

phytos.mat <- as.matrix(phytos)


# ---- NMDS ----


NMS <- metaMDS(phytos.mat, distance = "bray", k = 3)
goodness(NMS)
stressplot(NMS)

# plot(NMS, type = "t")

phytos.nmds <- as.data.frame(NMS$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
phytos.nmds$Date <- parse_date_time(rownames(phytos.nmds), orders = "ymd")  # create a column of site names, from the rownames of data.scores
phytos.nmds$Year <- year(phytos.nmds$Date)
phytos.nmds$Month <- month(phytos.nmds$Date)
phytos.nmds$Day <- day(phytos.nmds$Date)

# combo <- merge(phytos.nmds, aop.df, by = "Date", all.x = T, all.y = F)


phytos.species <- as.data.frame(NMS$species)
phytos.species$NMS.length <- sqrt(phytos.species$MDS1^2 + phytos.species$MDS2^2)
phytos.species.top <- phytos.species[head(order(phytos.species$NMS.length, decreasing = T), n = 100),]

# phytos.species.top$seq <- rownames(phytos.species.top)
# map.all <- rbind(map.arc, map.bac, map.euk)
# colnames(map.all)[1] <- "seq"
# phytos.species.top <- merge(x = phytos.species.top, y = map.all, by = "seq", all.x = T, all.y = F)
# 
# taxa.bac <- taxa.bac[,which(colnames(taxa.bac) %in% c("X", "taxon"))]
# taxa.arc <- taxa.arc[,which(colnames(taxa.arc) %in% c("X", "taxon"))]
# taxa.euk <- taxa.euk[,which(colnames(taxa.euk) %in% c("X", "taxon"))]
# taxa.all <- rbind(taxa.arc, taxa.bac, taxa.euk)
# 
# colnames(taxa.all)[1] <- "global_edge_num"
# phytos.species.top <- merge(phytos.species.top, y = taxa.all, by = "global_edge_num", all.x = T, all.y = F)
# 
# try.it <- phytos.species.top
# 
# 
# scale.factor <- 1.2

# AOP
ggplot() +
  geom_point(data = phytos.nmds, aes(x = MDS1, y = MDS2), alpha = 1, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  # geom_segment(data = phytos.species.top[head(order(phytos.species.top$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = MDS1/3, yend = MDS2/3), color = "black") +
  # geom_text(data = phytos.species.top[head(order(phytos.species.top$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(MDS1*1.2/3, factor = 100), y = jitter(MDS2*1.2/3, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phytos.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phytos.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phytos.species.model.predictors)), color = "green4") +
  labs(color = "AOP", size = "abs(AOP)", x = "Dim 1", y = "Dim 2") +
  ggtitle("NMDS: W. E. Allen Phytoplankton Microscopy Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 





