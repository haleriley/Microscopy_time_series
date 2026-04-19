# identifying differentially abundant taxa between time series
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
library(gplots)
library(viridis)


# ---- read in data ----

setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/")

# combo.full <- readRDS("2026-02-09_microscopy_o2bio_combo_full.rds")
# combo.full <- readRDS("2026-03-23_microscopy_o2bio_combo_full.rds")
combo.full <- readRDS("2026-03-27_microscopy_o2bio_combo_full.rds")
# combo.sio.temp <- readRDS("2026-02-14_combo_sio_temp.rds") ## heatwaves only outside climatological baseline
combo.sio.temp <- readRDS("2026-02-14_combo_sio_temp.rds") ## all heatwaves

model.var.imp <- readRDS("2026-02-09_microscopy_model_var_imp.rds")

# predictors <- readRDS("2026-02-09_microscopy_model_predictors.rds")
# predictors <- readRDS("2026-03-23_phyto_rf_predictors.rds")
predictors <- readRDS("2026-03-27_phyto_rf_predictors.rds")
# predictors <- predictors[which(substr(predictors, 1, 5) != "other")]

unique.dinos <- readRDS("2026-02-09_unique_dinos.rds")
unique.diatoms <- readRDS("2026-02-09_unique_diatomss.rds")

model.var.imp <- readRDS("2026-02-09_microscopy_model_var_imp.rds")

# # ---- quick time series of rel abund ----
# 
# combo.full.long <- combo.full %>% pivot_longer(cols = all_of(predictors), names_to = "Taxon", values_to = "Rel.Abund")
# ggplot(combo.full.long) +
#   # geom_line(aes(x = Date, y = Rel.Abund, color = Taxon)) +
#   geom_area(aes(x = Date, y = Rel.Abund, fill = Taxon)) +
#   theme_bw()

# ---- optional: hellinger transformation of phyto rel abunds ----

combo.full <- combo.full[,-which(colnames(combo.full) == "O2bio.estimated")]
combo.full <- na.omit(combo.full)
combo.full[,predictors] <- decostand(combo.full[,predictors], method = "hellinger")


# ---- boxplot comparison of abundances ----

# combo.full$Dataset[which(substr(combo.full$Dataset, 1, 1) == "W")] <- "WE Allen"
combo.full$Blob.status <- "Blob"
combo.full$Blob.status[which(year(combo.full$Date) <= 2013)] <- "PreBlob"
combo.full$Blob.status[which(year(combo.full$Date) >= 2016)] <- "PostBlob"

combo.full.long <- combo.full %>% pivot_longer(cols = all_of(predictors), names_to = "Taxon", values_to = "Count")

combo.full.long$Group <- NA
combo.full.long$Group[which(combo.full.long$Taxon %in% unique.diatoms)] <- "Diatom"
combo.full.long$Group[which(combo.full.long$Taxon %in% unique.dinos)] <- "Dinoflagellate"

combo.full.long$Blob.status <- factor(combo.full.long$Blob.status, levels = c("PreBlob", "Blob", "PostBlob"))
combo.full.long$Taxon <- factor(combo.full.long$Taxon)

ggplot(data = combo.full.long) +
  geom_boxplot(aes(x = Taxon, y = Count, fill = Blob.status)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(data = combo.full.long) +
  geom_boxplot(aes(x = Blob.status, y = Count, fill = Blob.status)) +
  facet_wrap(.~Taxon, scales = "free") +
  theme_bw()

combo.full.long.dinos <- combo.full.long[which(combo.full.long$Taxon %in% unique.dinos),]
combo.full.long.diatoms <- combo.full.long[which(combo.full.long$Taxon %in% unique.diatoms),]

ggplot(data = combo.full.long.dinos) +
  geom_boxplot(aes(x = Taxon, y = log10(Count), fill = Blob.status)) +
  theme_bw() +
  scale_fill_manual(values = c("blue", "red", "yellow")) +
  theme(axis.text.x = element_text(angle = 90))

ggplot(data = combo.full.long.diatoms) +
  geom_boxplot(aes(x = Taxon, y = log10(Count), fill = Blob.status)) +
  theme_bw() +
  scale_fill_manual(values = c("blue", "red", "yellow")) +
  theme(axis.text.x = element_text(angle = 90))

# ---- heatmap of phytos ----

combo.full.mat <- as.matrix(na.omit(combo.full[,which(colnames(combo.full) %in% predictors)]))
rownames(combo.full.mat) <- as.character(combo.full$Date)

heat.cols <- viridis::viridis(100, direction = 1)

sample.colors <- c(rep("red", times = length(which(substr(rownames(combo.full.mat), start = 1, stop = 1) == "1"))), 
                   rep("gold", times = length(which(substr(rownames(combo.full.mat), start = 1, stop = 1) == "2"))))

tax.colors <- rep(NA, times = length(predictors))
tax.colors[which(colnames(combo.full.mat) %in% unique.dinos)] <- "darkgreen"
tax.colors[which(colnames(combo.full.mat) %in% unique.diatoms)] <- "blue"

dev.off()
heatmap.2(t(combo.full.mat),
          trace = 'none',
          #Colv = NA,
          scale = NULL,
          col = heat.cols,
          # labRow = tally.lab.Row[selected],
          # margins = c(10,10),
          colCol = sample.colors,
          colRow = tax.colors,
          key = TRUE,
          keysize = 2,
          cexRow = 0.5,
          cexCol = 0.2,
)




# ---- identify differentially abundant taxa using DEseq2 ----

## lots of this code is from 2023-11-06_deseq(1).R


combo.full.mat <- combo.full[,which(colnames(combo.full) %in% predictors)]
combo.full.mat <- combo.full.mat

rownames(combo.full.mat) <- combo.full$Date


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

metadata <- data.frame(Date = combo.full$Date)
metadata$Blob.status <- "Blob"
metadata$Blob.status[which(year(metadata$Date) <= 2013)] <- "PreBlob"
metadata$Blob.status[which(year(metadata$Date) >= 2016)] <- "PostBlob"
my.rownames <- metadata$Date
metadata <- as.data.frame(metadata[,-1])
rownames(metadata) <- my.rownames
colnames(metadata) <- "Blob.status"

# metadata$Dataset[which(substr(metadata$Dataset, 1,1) == "W")] <- "WEAllen"
# metadata$heatwave.perc.90 <- factor(metadata$Dataset)
# rownames(metadata) <- combo.full$Date

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
library(ggplot2)
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


res_tax_sig <- res_tax_sig[order(res_tax_sig$log2FoldChange, decreasing = T),]
res_tax_sig$Predictors <- factor(res_tax_sig$OTU, levels = res_tax_sig$OTU[order(res_tax_sig$log2FoldChange, decreasing = T)])

ggplot(data = res_tax_sig) +
  geom_bar(aes(x = Predictors, y = log2FoldChange, fill = log2FoldChange), stat = "identity") +
  theme_bw() +
  scale_fill_viridis_c() +
  coord_flip()


# colnames(res_tax_sig)[which(colnames(res_tax_sig) == "OTU")] <- "Predictors"

taxa.df <- merge(model.var.imp, res_tax_sig, by = "Predictors", all = T)
taxa.df$Predictors <- factor(taxa.df$Predictors, levels = taxa.df$Predictors[order(taxa.df$log2FoldChange, decreasing = T)])

taxa.df$Group <- NA
taxa.df$Group[which(taxa.df$Predictors %in% unique.diatoms)] <- "Diatom"
taxa.df$Group[which(taxa.df$Predictors %in% unique.dinos)] <- "Dinoflagellate"

ggplot(data = taxa.df) +
  geom_bar(aes(x = Predictors, y = log2FoldChange, 
               # fill = Variable.Importance,
               fill = log2FoldChange
               ), stat = "identity") +
  theme_bw() +
  scale_fill_viridis_c() +
  coord_flip() +
  theme(axis.text.y = element_text(color = rev(model.var.imp$Group))) +
  annotate("segment", x = "bacteriastrum_spp", y = -3, xend = "bacteriastrum_spp", yend = -5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("segment", x = "bacteriastrum_spp", y = 3, xend = "bacteriastrum_spp", yend = 5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  # annotate("text", x = "dinophysis_acuminata", y = -4, label = "SCCOOS", size = 3) +
  # annotate("text", x = "dinophysis_acuminata", y = 4, label = "WEAllen", size = 3) +
  annotate("text", x = "chaetoceros_concavicornis", y = -4, label = "SCCOOS", size = 3) +
  annotate("text", x = "chaetoceros_concavicornis", y = 4, label = "WEAllen", size = 3) +
  # labs(fill = "Predictor \nModel Importance")
  labs(fill = "log2FoldChange")



# confirming which (+/-) is SCCOOS or WE Allen
# ggplot(data = combo.full) +
#   geom_boxplot(aes(x = Dataset, y = pseudonitzschia_delicatissima)) +
#   theme_bw()


# --- identify differentially abundant taxa using DESeq2 (heatwavees) ----


combo.full.mat <- combo.full[,which(colnames(combo.full) %in% predictors)]
combo.full.mat <- combo.full.mat*10000



#We will convert our table to DESeqDataSet object
countData = round(as(combo.full.mat, "matrix"), digits = 0)
countData[is.na(countData)] <- 0
# countData <- countData[,which(colSums(countData) >= 100)]

## We will add 1 to the countData otherwise DESeq will fail with the error:
## estimating size factors
## Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  :
## every gene contains at least one zero, cannot compute log geometric means
countData<-(t(countData+1)) 


metadata <- as.data.frame(combo.sio.temp$heatwave.perc.90[which(combo.sio.temp$Date %in% combo.full$Date)])
colnames(metadata) <- "heatwave.perc.90"

# metadata$Dataset[which(substr(metadata$Dataset, 1,1) == "W")] <- "WEAllen"
metadata$heatwave.perc.90 <- factor(metadata$heatwave.perc.90)


dds <- DESeqDataSetFromMatrix(countData, metadata, as.formula(~ heatwave.perc.90))


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
library(ggplot2)
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


res_tax_sig <- res_tax_sig[order(res_tax_sig$log2FoldChange, decreasing = T),]
res_tax_sig$OTU <- factor(res_tax_sig$OTU, levels = res_tax_sig$OTU[order(res_tax_sig$log2FoldChange, decreasing = T)])

ggplot(data = res_tax_sig) +
  geom_bar(aes(x = OTU, y = log2FoldChange, fill = log2FoldChange), stat = "identity") +
  theme_bw() +
  scale_fill_viridis_c() +
  coord_flip()


colnames(res_tax_sig)[which(colnames(res_tax_sig) == "OTU")] <- "Predictors"

taxa.df <- merge(model.var.imp, res_tax_sig, by = "Predictors")
taxa.df$Predictors <- factor(taxa.df$Predictors, levels = taxa.df$Predictors[order(taxa.df$log2FoldChange, decreasing = T)])

taxa.df$Group <- NA
taxa.df$Group[which(taxa.df$Predictors %in% unique.diatoms)] <- "Diatom"
taxa.df$Group[which(taxa.df$Predictors %in% unique.dinos)] <- "Dinoflagellate"

ggplot(data = taxa.df) +
  geom_bar(aes(x = Predictors, y = log2FoldChange, fill = log2FoldChange), stat = "identity") +
  theme_bw() +
  scale_fill_viridis_c() +
  coord_flip() +
  theme(axis.text.y = element_text(color = rev(model.var.imp$Group))) +
  annotate("segment", x = "bacteriastrum_spp", y = -3, xend = "bacteriastrum_spp", yend = -5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("segment", x = "bacteriastrum_spp", y = 3, xend = "bacteriastrum_spp", yend = 5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x = "guinardia_striata", y = -4, label = "No Heatwave", size = 3) +
  annotate("text", x = "guinardia_striata", y = 4, label = "Heatwave", size = 3) +
  theme(legend.position = "none")

# confirming which (+/-) is SCCOOS or WE Allen
ggplot(data = combo.sio.temp) +
  geom_boxplot(aes(x = heatwave.perc.90, y = ceratium_furca)) +
  theme_bw()



# --- identify differentially abundant taxa using DESeq2 (coldwavees) ----


combo.full.mat <- combo.full[,which(colnames(combo.full) %in% predictors)]
combo.full.mat <- combo.full.mat*10000



#We will convert our table to DESeqDataSet object
countData = round(as(combo.full.mat, "matrix"), digits = 0)
countData[is.na(countData)] <- 0
# countData <- countData[,which(colSums(countData) >= 100)]

## We will add 1 to the countData otherwise DESeq will fail with the error:
## estimating size factors
## Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  :
## every gene contains at least one zero, cannot compute log geometric means
countData<-(t(countData+1)) 


metadata <- as.data.frame(combo.sio.temp$coldwave.perc.90[which(combo.sio.temp$Date %in% combo.full$Date)])
colnames(metadata) <- "coldwave.perc.90"

# metadata$Dataset[which(substr(metadata$Dataset, 1,1) == "W")] <- "WEAllen"
metadata$coldwave.perc.90 <- factor(metadata$coldwave.perc.90)


dds <- DESeqDataSetFromMatrix(countData, metadata, as.formula(~ coldwave.perc.90))


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
library(ggplot2)
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


res_tax_sig <- res_tax_sig[order(res_tax_sig$log2FoldChange, decreasing = T),]
res_tax_sig$OTU <- factor(res_tax_sig$OTU, levels = res_tax_sig$OTU[order(res_tax_sig$log2FoldChange, decreasing = T)])

ggplot(data = res_tax_sig) +
  geom_bar(aes(x = OTU, y = log2FoldChange, fill = log2FoldChange), stat = "identity") +
  theme_bw() +
  scale_fill_viridis_c() +
  coord_flip()


colnames(res_tax_sig)[which(colnames(res_tax_sig) == "OTU")] <- "Predictors"

taxa.df <- merge(model.var.imp, res_tax_sig, by = "Predictors")
taxa.df$Predictors <- factor(taxa.df$Predictors, levels = taxa.df$Predictors[order(taxa.df$log2FoldChange, decreasing = T)])

taxa.df$Group <- NA
taxa.df$Group[which(taxa.df$Predictors %in% unique.diatoms)] <- "Diatom"
taxa.df$Group[which(taxa.df$Predictors %in% unique.dinos)] <- "Dinoflagellate"

ggplot(data = taxa.df) +
  geom_bar(aes(x = Predictors, y = log2FoldChange, fill = log2FoldChange), stat = "identity") +
  theme_bw() +
  scale_fill_viridis_c() +
  coord_flip() +
  theme(axis.text.y = element_text(color = rev(model.var.imp$Group))) +
  annotate("segment", x = "bacteriastrum_spp", y = -3, xend = "bacteriastrum_spp", yend = -5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("segment", x = "bacteriastrum_spp", y = 3, xend = "bacteriastrum_spp", yend = 5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x = "guinardia_striata", y = -4, label = "No coldwave", size = 3) +
  annotate("text", x = "guinardia_striata", y = 4, label = "coldwave", size = 3) +
  theme(legend.position = "none")

# confirming which (+/-) is SCCOOS or WE Allen
ggplot(data = combo.sio.temp) +
  geom_boxplot(aes(x = coldwave.perc.90, y = ceratium_furca)) +
  theme_bw()



# ---- look at phytoplankton count community structure through NMDS ----

com.structure.for.nmds <- combo.full[,which(colnames(combo.full) %in% c(predictors, "Date"))]

## note that these are only predictors with model importance > 2%
ggplot(data = model.var.imp[order(model.var.imp$Variable.Importance, decreasing = T),]) +
  geom_bar(aes(x = Predictors, y = Variable.Importance), stat= "identity") +
  theme_bw() +
  coord_flip()

my.top.pred.taxa <- as.character(model.var.imp$Predictors[which(model.var.imp$Variable.Importance >= 5)])

com.structure.for.nmds <- combo.full[,which(colnames(combo.full) %in% c(my.top.pred.taxa, "Date"))]


## for whatever reason, something is fucky with 1924-03-17
# com.structure.for.nmds <- com.structure.for.nmds[which(com.structure.for.nmds$Date != parse_date_time("1924-03-17", orders = "Ymd")),]


# rm.dups.index <- which(duplicated(com.structure.for.nmds$Date))
# com.structure.for.nmds <- com.structure.for.nmds[-c((rm.dups.index)-1),]
my.rownames <- com.structure.for.nmds$Date
com.structure.for.nmds <- com.structure.for.nmds[,-1]
rownames(com.structure.for.nmds) <- my.rownames

com.structure.for.nmds <- na.omit(com.structure.for.nmds)

com.structure.for.nmds <- com.structure.for.nmds[which(rowSums(com.structure.for.nmds) != 0),]


# #### TEMPORARILY REMOVE TAXA THAT DON'T APPEAR IN WEALLEN TIME SERIES ----
# com.structure.for.nmds$Dataset <- "W.E. Allen"
# com.structure.for.nmds$Dataset[which(parse_date_time(rownames(com.structure.for.nmds), orders = "Ymd") >= parse_date_time(2000, orders = "Y"))] <- "SCCOOS"
# 
# weallen.phytos <- com.structure.for.nmds[which(com.structure.for.nmds$Dataset == "W.E. Allen"),]
# sccoos.phytos <- com.structure.for.nmds[which(com.structure.for.nmds$Dataset == "SCCOOS"),]
# 
# no.weallen.phytos <- names(which(colSums(weallen.phytos[,-ncol(weallen.phytos)]) == 0))
# # which(colSums(sccoos.phytos[,c(9:84)]) == 0)
# 
# com.structure.for.nmds <- com.structure.for.nmds[,-which(colnames(com.structure.for.nmds) %in% no.weallen.phytos)]
# com.structure.for.nmds <- com.structure.for.nmds[,which(colnames(com.structure.for.nmds) != "Dataset")]
#### ----



phyto.mat <- as.matrix(com.structure.for.nmds)
# phyto.mat <- phyto.mat
# phyto.dist <- dist(phyto.mat)

# dimcheckMDS(phyto.mat)



NMS <- metaMDS(phyto.mat, distance = "bray", k = 4, autotransform = T)
goodness(NMS)
stressplot(NMS)

# plot(NMS)
# data.scores = as.data.frame(scores(NMS)$sites)
# plot(data.scores[,1:2])

phyto.nmds <- as.data.frame(NMS$points)


## reducing the number of taxa for envfit analysis (assuming only abundant taxa influence community structure strongly)
summary(colSums(phyto.mat))
abundant.taxa.index <- which(colSums(phyto.mat) >= 1)
phyto.mat <- phyto.mat[,abundant.taxa.index]

print("It's Happening!!! Running envfit for phyto NMDS species scores")

phyto.envfit <- envfit(NMS,phyto.mat,nperm=0)
saveRDS(phyto.envfit, file = "2026-03-27_phyto_envfit_nmds.rds")

phyto.envfit <- readRDS("2026-03-27_phyto_envfit_nmds.rds")

# try.it <- phyto.envfit$vectors
phyto.taxa.nmds <- as.data.frame(scores(phyto.envfit, "vectors"))
# phyto.taxa.nmds <- as.data.frame(phyto.envfit$vectors$arrows)
# phyto.taxa.nmds.r <- phyto.envfit$vectors$r

phyto.species.save <- phyto.taxa.nmds
saveRDS(phyto.species.save, file = "2026-03-27_phyto_species.rds")

# try.it <- as.data.frame(scores.envfit(phyto.envfit))
# saveRDS(try.it, "2025-11-04_try_it.rds")

# phyto.nmds <- as.data.frame(NMS$points)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
phyto.nmds$Date <- parse_date_time(rownames(phyto.nmds), orders = "ymd")  # create a column of site names, from the rownames of data.scores
phyto.nmds$Year <- year(phyto.nmds$Date)
phyto.nmds$Month <- month(phyto.nmds$Date)
phyto.nmds$Day <- day(phyto.nmds$Date)

combo <- merge(phyto.nmds, combo.full, by = "Date", all.x = T, all.y = F)

print("Combining was a success!")


phyto.species <- phyto.taxa.nmds
# phyto.species <- as.data.frame(scores(NMS)$species)
# plot(phyto.species[,1:2])
phyto.species$taxon <- rownames(phyto.species)
# phyto.species <- get.taxon.info(phyto.species)
phyto.species$NMS.length <- sqrt(phyto.species$NMDS1^2 + phyto.species$NMDS2^2)
# colnames(phyto.species)[3:4] <- c("MDS1", "MDS2")
summary(phyto.species$NMS.length)
# phyto.species.top <- phyto.species[which(phyto.species$NMS.length >= 0.4),]
# phyto.species.top <- phyto.species[head(order(phyto.species$NMS.length, decreasing = T), n = round(0.1*nrow(phyto.species))),]
# nrow(phyto.species.top)/ncol(phyto.species)

# saveRDS(phyto.species.top, file = "2026-03-21_phyto_species_top.rds")

# print("Saving species was a success!")

phyto.nmds$Dataset <- "W.E. Allen"
phyto.nmds$Dataset[which(parse_date_time(rownames(phyto.nmds), orders = "Ymd") >= parse_date_time(2000, orders = "Y"))] <- "SCCOOS"

com.structure.for.nmds$Dataset <- "W.E. Allen"
com.structure.for.nmds$Dataset[which(parse_date_time(rownames(com.structure.for.nmds), orders = "Ymd") >= parse_date_time(2000, orders = "Y"))] <- "SCCOOS"

# weallen.phytos <- com.structure.for.nmds[which(com.structure.for.nmds$Dataset == "W.E. Allen"),]
sccoos.phytos <- com.structure.for.nmds[which(com.structure.for.nmds$Dataset == "SCCOOS"),]

# no.weallen.phytos <- names(which(colSums(weallen.phytos[,-ncol(weallen.phytos)]) == 0))



# scale.factor <- 2


ggplot() +
  geom_point(data = phyto.nmds, aes(x = MDS1, y = MDS2, color = factor(Month)), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds, aes(x = MDS1, y = MDS2, color = factor(Month), shape = Dataset), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  # scale_color_viridis_d() +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  # stat_ellipse(data = phyto.nmds, aes(x = MDS1, y = MDS2, fill = factor(Month), color = factor(Month)), geom = "polygon", alpha = 0.1) +
  # scale_fill_viridis_d() +
  scale_color_manual(values = rainbow(12)) +
  scale_fill_manual(values = rainbow(12)) +
  geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Month", fill = "Month", x = "Dim 1", y = "Dim 2") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 

ggplot() +
  geom_point(data = phyto.nmds, aes(x = MDS1, y = MDS2, color = factor(Year)), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds, aes(x = MDS1, y = MDS2, color = factor(Month), shape = Dataset), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_d() +
  # geom_path(data = phyto.nmds, aes(x = MDS1, y = MDS2), size = 0.5, alpha = 0.3) +
  # stat_ellipse(data = phyto.nmds, aes(x = MDS1, y = MDS2, fill = factor(Year), color = factor(Year)), geom = "polygon", alpha = 0.2) +
  # scale_fill_viridis_d() +
  # scale_color_manual(values = rainbow(12)) +
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Year", fill = "Year", x = "Dim 1", y = "Dim 2") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 


# 
# for(y in 1:length(unique(phyto.nmds$Year))){
#   
#   my.year <- unique(phyto.nmds$Year)[y]
#   
#   my.df <- phyto.nmds[which(phyto.nmds$Year <= my.year),]
#   
#   png(file = paste("Figures/NMDS_phyto_shift/2026-04-03_nmds_", my.year, ".png", sep = ""), width = 700, height = 500)
#   
#   ggplot() +
#     geom_point(data = my.df, aes(x = MDS1, y = MDS2, color = factor(Year)), alpha = 1, size = 3) +
#     scale_color_manual(n = length(unique(phyto.nmds$Year))) +
#     stat_ellipse(data = my.df[which(my.df$Year == my.year),], aes(x = MDS1, y = MDS2, fill = factor(Year), color = factor(Year)), geom = "polygon", alpha = 0.2) +
#     scale_color_manual(viridis(n = length(unique(phyto.nmds$Year)))) +
#     # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
#     # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
#     labs(color = "Year", fill = "Year", x = "Dim 1", y = "Dim 2") +
#     ggtitle(paste(my.year)) +
#     theme_bw() +
#     theme(panel.grid = element_blank()) +
#     theme(axis.title = element_text(size = 14, face = "bold"), 
#           axis.text = element_text(size = 12), 
#           legend.text = element_text(size = 12), 
#           legend.title = element_text(size = 14, face = "bold"),
#           title = element_text(face = "bold")) +
#     ylim(c(1.1*min(phyto.nmds$MDS2), 1.1*max(phyto.nmds$MDS2))) +
#     xlim(c(1.1*min(phyto.nmds$MDS1), 1.1*max(phyto.nmds$MDS1))) 
#   
#   
#   dev.off()
#   
#     
#   
# }


com.structure.for.nmds.to.merge <- com.structure.for.nmds
com.structure.for.nmds.to.merge$Date <- parse_date_time(rownames(com.structure.for.nmds.to.merge), orders = "Ymd")

combo.taxa.dim1.comparison <- merge(phyto.nmds, com.structure.for.nmds.to.merge, by = "Date")
combo.taxa.dim2.comparison <- merge(phyto.nmds, com.structure.for.nmds.to.merge, by = "Date")


my.df <- data.frame("Taxon" = predictors, "Dim1.cor" = NA)

for(t in 1:length(colnames(sccoos.phytos))){
  
  my.taxon <- colnames(sccoos.phytos)[t]
  
  my.cor <- abs(cor(x = combo.taxa.dim1.comparison$MDS1, combo.taxa.dim1.comparison[,which(colnames(combo.taxa.dim1.comparison) == my.taxon)]))
  
  my.df[which(my.df$Taxon == my.taxon),2] <- my.cor
  
}

my.df <- data.frame("Taxon" = predictors, "Dim2.cor" = NA)

for(t in 1:length(colnames(sccoos.phytos))){
  
  my.taxon <- colnames(sccoos.phytos)[t]
  
  my.cor <- abs(cor(x = combo.taxa.dim1.comparison$MDS2, combo.taxa.dim1.comparison[,which(colnames(combo.taxa.dim1.comparison) == my.taxon)]))
  
  my.df[which(my.df$Taxon == my.taxon),2] <- my.cor
  
}

# ---- combine with heatwave data ----

# combo.sio.temp <- readRDS("2026-02-14_combo_sio_temp.rds") ## heatwaves only outside climatological baseline
combo.sio.temp <- readRDS("2026-04-14_combo_sio_temp.rds") ## all heatwaves
combo.sio.temp.to.merge <- combo.sio.temp[,c("SURF_TEMP_C", "heatwave.perc.90", "coldwave.perc.90", "Date")]

phyto.nmds.temp <- merge(x = phyto.nmds, y = combo.sio.temp.to.merge, by = "Date", all.x = T, all.y = F)


## colored by surface water temp
ggplot() +
  geom_point(data = phyto.nmds.temp, aes(x = MDS4, y = MDS2, shape = Dataset, color = SURF_TEMP_C), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_c() +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Dataset", x = "Dim 4", y = "Dim 2") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 


## colored by heatwaves
ggplot() +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "no"),], aes(x = MDS1, y = MDS2, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "yes"),], aes(x = MDS1, y = MDS2, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_d() +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Heatwave?", x = "Dim 1", y = "Dim 2") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 

ggplot() +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "no"),], aes(x = MDS1, y = MDS3, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "yes"),], aes(x = MDS1, y = MDS3, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_d() +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Heatwave?", x = "Dim 1", y = "Dim 3") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 


ggplot() +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "no"),], aes(x = MDS1, y = MDS4, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "yes"),], aes(x = MDS1, y = MDS4, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_d() +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Heatwave?", x = "Dim 1", y = "Dim 4") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 

ggplot() +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "no"),], aes(x = MDS1, y = MDS4, color = factor(Year)), alpha = 1, size = 1) +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "yes"),], aes(x = MDS1, y = MDS4, color = factor(Year)), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_d() +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Heatwave?", x = "Dim 1", y = "Dim 4") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 

summary(lm(phyto.nmds.temp$SURF_TEMP_C~phyto.nmds.temp$MDS1))
summary(lm(phyto.nmds.temp$SURF_TEMP_C~phyto.nmds.temp$MDS4))
summary(lm(phyto.nmds.temp$Year~phyto.nmds.temp$MDS1))
summary(lm(phyto.nmds.temp$Year~phyto.nmds.temp$MDS4))

phyto.nmds.temp$Blob.status <- "Blob"
phyto.nmds.temp$Blob.status[which(phyto.nmds.temp$Year <= 2013)] <- "PreBlob"
phyto.nmds.temp$Blob.status[which(phyto.nmds.temp$Year >= 2016)] <- "PostBlob"


ggplot() +
  geom_point(data = phyto.nmds.temp, aes(x = MDS1, y = MDS2, color = Blob.status, size = heatwave.perc.90), alpha = 1) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_d() +
  scale_size_manual(values = c(2,4)) +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Blob Status", size = "Heatwave?", x = "Dim 1", y = "Dim 2") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 




dev.off()
for(y in 1:length(unique(phyto.nmds.temp$Year))){
  
  my.year <- unique(phyto.nmds.temp$Year)[y]
  
  my.df <- phyto.nmds.temp[which(phyto.nmds.temp$Year <= my.year),]
  
  my.colors <- viridis(n = length(unique(phyto.nmds.temp$Year)))
  
  png(file = paste("Figures/NMDS_phyto_shift/2026-04-06_nmds_year_Dim1+2_nootherdiatomsdinos", my.year, ".png", sep = ""), width = 700, height = 500)
  
  a <- ggplot() +
    geom_point(data = my.df, aes(x = MDS1, y = MDS2, color = factor(Year)), alpha = 1, size = 3) +
    # scale_color_manual(n = length(unique(phyto.nmds$Year))) +
    stat_ellipse(data = my.df[which(my.df$Year == my.year),], aes(x = MDS1, y = MDS2, fill = factor(Year), color = factor(Year)), geom = "polygon", alpha = 0.2) +
    scale_color_manual(values = my.colors[1:y]) +
    scale_fill_manual(values = my.colors[y]) +
    # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
    # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
    labs(color = "Year", fill = "Year", x = "Dim 1", y = "Dim 2") +
    ggtitle(paste(my.year)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(axis.title = element_text(size = 14, face = "bold"), 
          axis.text = element_text(size = 12), 
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 14, face = "bold"),
          title = element_text(face = "bold")) +
    ylim(c(1.1*min(phyto.nmds$MDS2), 1.1*max(phyto.nmds$MDS2))) +
    xlim(c(1.1*min(phyto.nmds$MDS1), 1.1*max(phyto.nmds$MDS1))) 
  
  print(a)
  
  dev.off()
  
  
  my.colors <- viridis(n = length(unique(phyto.nmds.temp$Blob.status)))
  my.colors.blob <- c(rep(my.colors[1], 3), rep(my.colors[2], 2), rep(my.colors[3], 9))
  
  png(file = paste("Figures/NMDS_phyto_shift/2026-04-06_nmds_blob_Dim1+2_nootherdiatomsdinos", my.year, ".png", sep = ""), width = 700, height = 500)
  
  b <- ggplot() +
    geom_point(data = my.df, aes(x = MDS1, y = MDS2, color = factor(Year)), alpha = 1, size = 3) +
    # scale_color_manual(n = length(unique(phyto.nmds$Year))) +
    stat_ellipse(data = my.df[which(my.df$Year == my.year),], aes(x = MDS1, y = MDS2, fill = factor(Year), color = factor(Year)), geom = "polygon", alpha = 0.2) +
    scale_color_manual(values = my.colors.blob[1:y]) +
    scale_fill_manual(values = my.colors.blob[y]) +
    # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
    # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
    labs(color = "Year", fill = "Year", x = "Dim 1", y = "Dim 2") +
    ggtitle(paste(my.year)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(axis.title = element_text(size = 14, face = "bold"), 
          axis.text = element_text(size = 12), 
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 14, face = "bold"),
          title = element_text(face = "bold")) +
    ylim(c(1.1*min(phyto.nmds$MDS2), 1.1*max(phyto.nmds$MDS2))) +
    xlim(c(1.1*min(phyto.nmds$MDS1), 1.1*max(phyto.nmds$MDS1))) +
    guides(fill = "none")
  
  print(b)
  
  dev.off()
  
  
}


heat.cols <- viridis::viridis(100, direction = 1)

phyto.mat.log10 <- log10(phyto.mat)
phyto.mat.log10[which(phyto.mat.log10 < 0)] <- 0
phyto.mat.log10 <- t(phyto.mat.log10)

my.colors.blob.2022 <- viridis(4)
my.colCol <- rep(my.colors.blob.2022[3], times = ncol(phyto.mat.log10))
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) <= 2013)] <- my.colors.blob.2022[1]
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) %in% c(2014, 2015))] <- my.colors.blob.2022[2]
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) >= 2022)] <- my.colors.blob.2022[4]


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
          adjCol = 0.5,
          Colv = F,
          colCol = my.colCol
          # key = TRUE,
          # keysize = 2,
          # density.info = "density"
)




ggplot() +
  geom_point(data = phyto.nmds.temp, aes(x = MDS1, y = MDS4, color = Blob.status, size = heatwave.perc.90), alpha = 1) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_d() +
  scale_size_manual(values = c(2,4)) +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Blob Status", size = "Heatwave?", x = "Dim 1", y = "Dim 4") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 



ggplot() +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "no"),], aes(x = MDS3, y = MDS2, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "yes"),], aes(x = MDS3, y = MDS2, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_d() +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Heatwave?", x = "Dim 3", y = "Dim 2") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 



ggplot() +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "no"),], aes(x = MDS4, y = MDS2, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "yes"),], aes(x = MDS4, y = MDS2, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_d() +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Heatwave?", x = "Dim 4", y = "Dim 2") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 

ggplot() +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "no"),], aes(x = MDS3, y = MDS4, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  geom_point(data = phyto.nmds.temp[which(phyto.nmds.temp$heatwave.perc.90 == "yes"),], aes(x = MDS3, y = MDS4, shape = Dataset, color = heatwave.perc.90), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_d() +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Heatwave?", x = "Dim 3", y = "Dim 4") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 


# phyto.nmds.temp.cor <- na.omit(phyto.nmds.temp[,c("MDS1", "MDS2", "MDS3", "MDS4", "SURF_TEMP_C")])
summary(lm(phyto.nmds.temp$SURF_TEMP_C~phyto.nmds.temp$MDS1))
summary(lm(phyto.nmds.temp$SURF_TEMP_C~phyto.nmds.temp$MDS2))
summary(lm(phyto.nmds.temp$SURF_TEMP_C~phyto.nmds.temp$MDS3))
summary(lm(phyto.nmds.temp$SURF_TEMP_C~phyto.nmds.temp$MDS4))




# ---- plot time series of most diff abundant MDS1 ----

phytos.combo <- merge(phyto.nmds.temp, com.structure.for.nmds.to.merge, by = c("Date", "Dataset"))


my.Dim1.df <- data.frame(taxon = rep(NA, length(predictors)), adj.R2 = rep(NA, length(predictors)), p.val = rep(NA, length(predictors)))

for(t in 1:length(predictors)){
  
  my.taxon <- predictors[t]
  
  temp <- summary(lm(phytos.combo[,which(colnames(phytos.combo) == my.taxon)]~phytos.combo$MDS1))
  
  my.Dim1.df$taxon[t] <- my.taxon
  my.Dim1.df$adj.R2[t] <- temp$adj.r.squared
  my.Dim1.df$p.val[t] <- temp$coefficients[2,4]
  
}

my.Dim1.df <- my.Dim1.df[which(my.Dim1.df$adj.R2 > 0 & my.Dim1.df$p.val < 0.05),]
my.Dim1.df <- my.Dim1.df[order(my.Dim1.df$adj.R2),]
my.Dim1.df$taxon <- factor(my.Dim1.df$taxon, levels = my.Dim1.df$taxon[order(my.Dim1.df$adj.R2)])

ggplot(data = my.Dim1.df[order(my.Dim1.df$adj.R2),]) +
  geom_bar(aes(x = taxon, y = adj.R2), stat = "identity") +
  theme_bw() +
  coord_flip()


phytos.combo <- merge(phytos.combo, combo[,c(1,124,125)], by = "Date")

summary(lm(phytos.combo$O2bio.predicted.altered~phytos.combo$MDS1))
summary(lm(phytos.combo$O2bio.predicted.altered~phytos.combo$MDS2))
summary(lm(phytos.combo$O2bio.predicted.altered~phytos.combo$MDS3))
summary(lm(phytos.combo$O2bio.predicted.altered~phytos.combo$MDS4))

ggplot() +
  geom_point(data = phytos.combo, aes(x = MDS1, y = MDS2, color = O2bio.predicted.altered), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_c() +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Dataset", x = "Dim 4", y = "Dim 2") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 

ggplot() +
  geom_point(data = phytos.combo, aes(x = MDS1, y = MDS4, color = O2bio.predicted.altered), alpha = 1, size = 3) +
  # geom_point(data = phyto.nmds[which(phyto.nmds$Year >= 2014),], aes(x = MDS1, y = MDS2), color = "red", alpha = 1, size = 3, shape = "+") +
  # geom_point(data = combo[which(is.na(combo$O2bio.estimated) == T),], aes(x = NMDS1, y = NMDS2), color = "grey69", alpha = 0.2, size = 3) +
  # scale_color_viridis_c() +
  # geom_text(aes(x = NMDS1, y = NMDS2, label = sample.date)) +
  # scale_color_manual(values = my.year.colors) +
  # scale_color_gradient2(low = "blue", high = "red", mid = "white") +
  scale_color_viridis_c() +
  # geom_path(aes(x = NMDS1, y = NMDS2), size = 0.5, alpha = 0.3) +
  
  # geom_segment(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), color = "black") +
  
  # geom_text(data = phyto.species[head(order(phyto.species$NMS.length, decreasing = TRUE), n = 10),], aes(x = jitter(NMDS1*1.2, factor = 100), y = jitter(NMDS2*1.2, factor = 200), label = taxon), color = "black") +
  # geom_segment(data = phyto.species.model.predictors, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), color = "green4") +
  # geom_text(data = phyto.species.model.predictors, aes(x = MDS1*scale.factor, y = MDS2*scale.factor, label = rownames(phyto.species.model.predictors)), color = "green4") +
  labs(color = "Predicted O2bio", x = "Dim 1", y = "Dim 2") +
  # ggtitle("NMDS: Microbial Community Time Series") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"),
        title = element_text(face = "bold")) 



# # ---- plot time series of most diff abundant and NMDS-influential taxa ----
# 
# my.top.taxa <- my.df$Taxon[which(my.df$Dim2.cor >= 0.3)]
# 
# ggplot(data = sccoos.combo.long[which(sccoos.combo.long$Taxa %in% my.top.taxa),]) +
#   geom_line(aes(x = Date, y = Count, color = Taxa)) +
#   facet_wrap(.~Taxa, ncol = 2, scales = "free") +
#   theme_bw()
# 

# ---- look at three shifts/four periods throughout time series - shifting taxa, NMDS Dim1, and 02bio values ----

## copied from above to remind:

heat.cols <- viridis::viridis(100, direction = 1)

phyto.mat.log10 <- log10(phyto.mat)
phyto.mat.log10[which(phyto.mat.log10 < 0)] <- 0
phyto.mat.log10 <- t(phyto.mat.log10)

summary(rowSums(phyto.mat.log10))
phyto.mat.log10 <- phyto.mat.log10[which(rowSums(phyto.mat.log10) >= 200),]
my.abund.taxa <- rownames(phyto.mat.log10)

my.colors.blob.2022 <- palette.colors(palette = "R4", n = 4)
my.colCol <- rep(my.colors.blob.2022[3], times = ncol(phyto.mat.log10))
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) <= 2013)] <- my.colors.blob.2022[1]
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) %in% c(2014, 2015))] <- my.colors.blob.2022[2]
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) >= 2022)] <- my.colors.blob.2022[4]


dev.off()
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
          adjCol = 0.5,
          Colv = F,
          colCol = my.colCol
          # key = TRUE,
          # keysize = 2,
          # density.info = "density"
)

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
          adjCol = 0.5,
          Colv = T,
          colCol = my.colCol
          # key = TRUE,
          # keysize = 2,
          # density.info = "density"
)


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
  theme_bw()

ggplot(phyto.df) +
  geom_line(aes(x = Date, y = shannon.div)) +
  geom_smooth(aes(x = Date, y = shannon.div), method = "loess") +
  theme_bw()

phyto.df$Blob.status <- "Blob"
phyto.df$Blob.status[which(year(phyto.df$Date) <= 2013)] <- "PreBlob"
phyto.df$Blob.status[which(year(phyto.df$Date) >= 2016)] <- "PostBlob"
phyto.df$Blob.status[which(year(phyto.df$Date) >= 2022)] <- "Post2022Shift"
phyto.df$Blob.status <- factor(phyto.df$Blob.status, levels = c("PreBlob", "Blob", "PostBlob", "Post2022Shift"))

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


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = richness, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Phytoplankton Richness") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
    inherit.aes = FALSE) +
  theme_bw()


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


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = shannon.div, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Phytoplankton Shannon Diversity") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
    inherit.aes = FALSE) +
  theme_bw()


## O2bio ##

phyto.df <- merge(phyto.df, combo[,which(colnames(combo) %in% c("Date", "O2bio.predicted.altered"))])

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


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = O2bio.predicted.altered, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Predicted [O2]bio") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
    inherit.aes = FALSE) +
  theme_bw()



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


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = SURF_TEMP_C, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Surface Water Temperature (C)") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
    inherit.aes = FALSE) +
  theme_bw()



# ---- measure volatility ----

sccoos <- read.csv("https://erddap.sccoos.org/erddap/tabledap/HABs-ScrippsPier.csv?time%2CChl_Volume_Filtered%2CChl1%2CChl2%2CAvg_Chloro%2CPhosphate%2CSilicate%2CNitrite%2CNitrite_Nitrate%2CAmmonium%2CNitrate%2CVolume_Settled_for_Counting&time%3E=2008-06-30T15%3A00%3A00Z&time%3C=2026-04-06T18%3A16%3A00Z")

sccoos <- sccoos[-1,]
colnames(sccoos)[1] <- c("Date")
sccoos$Date <- parse_date_time(sccoos$Date, orders= "Ymd HMS")
sccoos$Date <- parse_date_time(paste(year(sccoos$Date), month(sccoos$Date), day(sccoos$Date), sep = "-"), orders = "Ymd")
sccoos <- sccoos %>% mutate(across(c(colnames(sccoos)[c(2:ncol(sccoos))]), as.numeric))
sccoos <- sccoos %>% group_by(Date) %>% summarize_all(mean)
phyto.df <- merge(phyto.df, sccoos, by = "Date", all.x = T, all.y = F)


sccoos.FCM <- readRDS("../SCCOOS_microbial_time_series/R_Data/2026-03-17_sccoos_com_df_fcm_cleaned_rarefied_2000.rds")
sccoos.FCM <- sccoos.FCM[,which(colnames(sccoos.FCM) %in% c("Date", "total_ml_vol_corrected_AF", "total_ml_vol_corrected_SG"))]
phyto.df <- merge(phyto.df, sccoos.FCM, by = "Date", all.x = T, all.y = F)

phyto.df <- merge(phyto.df, combo.full[which(colnames(combo.full) %in% c("Date", "O2bio.estimated"))], all.x = T, all.y = F)

phyto.df$sum.phyto.counts <- rowSums(phyto.df[,which(colnames(phyto.df) %in% predictors)])



phyto.df$volatility.O2bio.predicted.altered <- NA
phyto.df$volatility.chlorophyll <- NA
phyto.df$volatility.fcm.AF <- NA
phyto.df$volatility.log10chlorophyll <- NA
phyto.df$volatility.O2bio.estimated <- NA
phyto.df$volatility.sum.count <- NA
phyto.df$volatility.log10.sum.count <- NA

for(d in 1:nrow(phyto.df)){
  
  my.date <- phyto.df$Date[d]
  
  my.df <- phyto.df[which(phyto.df$Date >= (my.date - 60*60*24*15) & phyto.df$Date <= (my.date + 60*60*24*15)),]
  
  phyto.df$volatility.O2bio.predicted.altered[d] <- sd(my.df$O2bio.predicted.altered, na.rm = T)
  phyto.df$volatility.chlorophyll[d] <- sd(my.df$Avg_Chloro, na.rm = T)
  phyto.df$volatility.log10chlorophyll[d] <- sd(log10(my.df$Avg_Chloro), na.rm = T)
  phyto.df$volatility.fcm.AF[d] <- sd(my.df$total_ml_vol_corrected_AF, na.rm = T)
  phyto.df$volatility.O2bio.estimated[d] <- sd(my.df$O2bio.estimated, na.rm = T)
  phyto.df$volatility.sum.count[d] <- sd(my.df$sum.phyto.count, na.rm = T)
  phyto.df$volatility.log10.sum.count[d] <- sd(log10(my.df$sum.phyto.counts), na.rm = T)
  
  
}

plot.vol.O2bio.predalt <- ggplot(data = phyto.df) +
  geom_line(aes(x = Date, y = volatility.O2bio.predicted.altered)) +
  theme_bw()

plot.vol.O2bio.est <- ggplot(data = phyto.df) +
  geom_line(aes(x = Date, y = volatility.O2bio.estimated)) +
  theme_bw()

# plot.vol.chlorophyll <- ggplot(data = phyto.df) +
#   geom_line(aes(x = Date, y = volatility.chlorophyll)) +
#   theme_bw()

plot.vol.log10chlorophyll <- ggplot(data = phyto.df) +
  geom_line(aes(x = Date, y = volatility.log10chlorophyll)) +
  theme_bw()

# plot.vol.FCM <- ggplot(data = phyto.df) +
#   geom_line(aes(x = Date, y = volatility.fcm.AF)) +
#   theme_bw()

# plot.vol.sum.counts <- ggplot(data = phyto.df) +
#   geom_line(aes(x = Date, y = volatility.sum.count)) +
#   theme_bw()

plot.vol.log10.sum.counts <- ggplot(data = phyto.df) +
  geom_line(aes(x = Date, y = volatility.log10.sum.count)) +
  theme_bw()

plot.vol.O2bio.est + plot.vol.log10chlorophyll + plot.vol.log10.sum.counts + plot_layout(ncol = 1)


summary(lm(phyto.df$volatility.log10chlorophyll~phyto.df$volatility.log10.sum.count))
ggplot(data = phyto.df) +
  geom_point(aes(x = volatility.log10.sum.count, y = volatility.log10chlorophyll)) +
  geom_smooth(aes(x = volatility.log10.sum.count, y = volatility.log10chlorophyll), method = "lm", se = F) +
  theme_bw()


## O2bio predicted ##
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
    aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
    inherit.aes = FALSE) +
  theme_bw()



## Chlorophyll ##
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
    aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
    inherit.aes = FALSE) +
  theme_bw()


## FCM counts ##
my.aov <- aov(phyto.df$volatility.fcm.AF~phyto.df$Blob.status)

my.tukey <- TukeyHSD(my.aov)

my.tukey.results <- my.tukey$`phyto.df$Blob.status`

my.df <- as.data.frame(my.tukey.results[which(my.tukey.results[,4] < 0.05), 4])
my.df$period.comparison <- rownames(my.df)

library(multcompView)

my.letters <- multcompLetters4(my.aov, my.tukey)[[1]]$Letters

y_max <- tapply(phyto.df$volatility.fcm.AF, phyto.df$Blob.status, max, na.rm = TRUE)

label_df <- data.frame(
  Blob.status = names(my.letters),
  Letter = my.letters,
  y_pos = y_max[names(my.letters)] + 0.2  # adjust spacing as needed
)


ggplot(phyto.df) +
  geom_boxplot(aes(x = Blob.status, y = volatility.fcm.AF, fill = Blob.status)) +
  scale_fill_manual(values = my.colors.blob.2022) +
  labs(x = "Time Period", y = "Volatility FCM Counts (AF)") +
  guides(fill = "none") +
  geom_text(
    data = label_df,
    aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
    inherit.aes = FALSE) +
  theme_bw()


summary(lm(phyto.df$volatility.O2bio.predicted.altered~phyto.df$volatility.log10chlorophyll))
ggplot(phyto.df) +
  geom_point(aes(x = volatility.log10chlorophyll, y = volatility.O2bio.predicted.altered)) +
  geom_smooth(aes(x = volatility.log10chlorophyll, y = volatility.O2bio.predicted.altered), method = "lm", se = F) +
  labs(x = "Volatility log10 Chlorophyll", y = "Volatility Predicted [O2]bio") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw()

summary(lm(phyto.df$volatility.O2bio.estimated~phyto.df$volatility.log10chlorophyll))
ggplot(phyto.df) +
  geom_point(aes(x = volatility.log10chlorophyll, y = volatility.O2bio.estimated)) +
  geom_smooth(aes(x = volatility.log10chlorophyll, y = volatility.O2bio.estimated), method = "lm", se = F) +
  labs(x = "Volatility log10 Chlorophyll", y = "Volatility Estimated [O2]bio - Model 1") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw()

summary(lm(phyto.df$volatility.O2bio.estimated~phyto.df$volatility.fcm.AF))
ggplot(phyto.df) +
  geom_point(aes(x = volatility.fcm.AF, y = volatility.O2bio.estimated)) +
  geom_smooth(aes(x = volatility.fcm.AF, y = volatility.O2bio.estimated), method = "lm", se = F) +
  labs(x = "Volatility log10 Chlorophyll", y = "Volatility Estimated [O2]bio - Model 1") +
  guides(fill = "none") +
  # geom_text(
  #   data = label_df,
  #   aes(x = Blob.status, y = log10(y_pos*2), label = Letter),
  #   inherit.aes = FALSE) +
  theme_bw()


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

summary(lm(phyto.df$volatility.log10chlorophyll~phyto.df$SURF_TEMP_C))
ggplot(data = phyto.df) +
  geom_point(aes(x = SURF_TEMP_C, y = volatility.log10chlorophyll)) +
  geom_smooth(aes(x = SURF_TEMP_C, y = volatility.log10chlorophyll), method = "lm", se = F) +
  theme_bw()


# ---- taxa variability between time periods ----# ---- tSURF_TEMP_Caxa variability between time periods ----

phytos.combo$Blob.status[which(phytos.combo$Year >= 2022)] <- "Post2022Shift"
phytos.combo$Blob.status <- factor(phytos.combo$Blob.status, levels = c("PreBlob", "Blob", "PostBlob", "Post2022Shift"))

library(ggsignif)
library(multcompView)

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

metadata <- data.frame(Date = combo.full$Date)
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
library(ggplot2)
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

ggplot(data = res_tax_sig) +
  geom_bar(aes(x = Predictors, y = log2FoldChange, fill = log2FoldChange), stat = "identity") +
  theme_bw() +
  scale_fill_viridis_c() +
  coord_flip()

# try.it <- phytos.combo[,which(colnames(phytos.combo) %in% c("lingulodinium_polyedra", "Blob.status"))]
# ggplot(try.it) +
  # geom_boxplot(aes(x = Blob.status, y = log10(lingulodinium_polyedra))) +
  # theme_bw()



## copied from above to remind:

heat.cols <- viridis::viridis(100, direction = 1)

phyto.mat.log10 <- log10(phyto.mat)
phyto.mat.log10[which(phyto.mat.log10 < 0)] <- 0
phyto.mat.log10 <- t(phyto.mat.log10)

summary(rowSums(phyto.mat.log10))
phyto.mat.log10 <- phyto.mat.log10[which(rowSums(phyto.mat.log10) >= 200),]
my.abund.taxa <- rownames(phyto.mat.log10)

my.colors.blob.2022 <- palette.colors(palette = "R4", n = 4)
my.colCol <- rep(my.colors.blob.2022[3], times = ncol(phyto.mat.log10))
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) <= 2013)] <- my.colors.blob.2022[1]
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) %in% c(2014, 2015))] <- my.colors.blob.2022[2]
my.colCol[which(as.numeric(substr(colnames(phyto.mat.log10), start = 1, stop = 4)) >= 2022)] <- my.colors.blob.2022[4]

my.rowCol <- rep("white", nrow(phyto.mat.log10))
increase.with.shift.taxa <- res_tax_sig$Predictors[which(res_tax_sig$log2FoldChange >= 2)]
my.rowCol[which(rownames(phyto.mat.log10) %in% increase.with.shift.taxa)] <- my.colors.blob.2022[4]   
decrease.with.shift.taxa <- res_tax_sig$Predictors[which(res_tax_sig$log2FoldChange <= -2)]
my.rowCol[which(rownames(phyto.mat.log10) %in% decrease.with.shift.taxa)] <- my.colors.blob.2022[3]   

heatmap.2(phyto.mat.log10,
          trace = 'none',
          dendrogram = "row",
          scale = NULL,
          col = heat.cols,
          cexRow = 0.7,
          cexCol = 0.5,
          key = FALSE,
          margins = c(3, 8),
          srtCol = 90,
          lhei = c(1,100),
          offsetCol = 0,
          adjCol = 0.5,
          # adjRow = 0.5,
          Colv = F,
          colCol = my.colCol,
          colRow = my.rowCol
          # key = TRUE,
          # keysize = 2,
          # density.info = "density"
)
  







