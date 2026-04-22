# Ecological Trajectory Analysis
# Original code from Isi
# RJH
# 2026-04-18



# ---- load libraries ----
library(tidyverse)
library(ranger)
library(lubridate)
library(vegan)
library(phyloseq)
library(scales)
library(grid)
library(reshape2)
library(compositions)
library(ggpubr)
library(tidyverse)
library(ecotraj)
library(tibble)

setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/R_Data")


combo.full <- readRDS("2026-03-27_microscopy_o2bio_combo_full.rds")
combo.sio.temp <- readRDS("2026-02-14_combo_sio_temp.rds") ## all heatwaves
model.var.imp <- readRDS("2026-02-09_microscopy_model_var_imp.rds")
predictors <- readRDS("2026-03-27_phyto_rf_predictors.rds")


setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/")



r_scripps <- combo.full[,which(colnames(combo.full) %in% c("Date", predictors))]
rownames(r_scripps) <- r_scripps$Date

## REMOVE 2021 BECAUSE IT'S MISSING A HUGE CHUNK OF DATA ##
r_scripps <- r_scripps[which(year(r_scripps$Date) != 2021),]

# r_scripps <- r_scripps[,-1]

# ##### transform to relative abundance and to Hellinger 
# #####https://github.com/joey711/phyloseq/issues/585
# ###### Normalize transforming to relative abundance from rarefied even depth dataset (abundance counts to fractional abundance as done in https://joey711.github.io/phyloseq/preprocess.html)
# head(sample_sums(r_scripps))
# r_hellinger_scripps <- transform_sample_counts(r_scripps, function(x) sqrt(x / sum(x)))

## skipping this step, just using absolute count data here, renaming this to not have to change the code too much
r_hellinger_scripps <- r_scripps


###################################################
#################SCRIPPS 2011 to 2024 ETA (ALL)
###################################################


# make a data frame from the sample_data
samples_df_scripps <- data.frame(sample_data(r_hellinger_scripps))

#assigning year and month
samples_df_scripps$day<-yday(samples_df_scripps$Date)  ###assigning day of the year

samples_df_scripps$week<-week(samples_df_scripps$Date)  ###assigning week day of the year

samples_df_scripps$year<-year(samples_df_scripps$Date)


metadata_scripps <- data.frame(year = year(samples_df_scripps$Date), month = month(samples_df_scripps$Date), day = yday(samples_df_scripps$Date))
rownames(metadata_scripps) <- rownames(samples_df_scripps)



######checking number of samples collected over time


#####Scripps
samples_df_scripps$point<- "1"

# table(samples_df_scripps$year)
# samples_df_scripps %>% count(samplesyear)


samples_scripps_all <- ggplot(samples_df_scripps, aes(x = day, y = point)) +
  facet_grid(rows = vars(year), scales = "free", space= "free_x") +
  geom_point(color = "#009E73") +
  theme(axis.text.x = element_text(face="plain", color="black", 
                                   size=5, angle=45, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", 
                                   size=10, angle=0),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"))+
  ylab("sample collection") +
  ggtitle("Scripps Pier time series") 
samples_scripps_all

ggsave(file="samples_scripps_all2018-2024.png", samples_scripps_all, path = "Figures/", width = 8, height = 5, units = "in")


metadata_scripps <- samples_df_scripps[, c("year","day")]



# Let us first define the vectors that describe the ecological entity and the survey of each observation:

names(metadata_scripps) <- c("entities", "surveys")

entities_scripps <- as.character(metadata_scripps$entities)
surveys_scripps <- as.integer(metadata_scripps$surveys)
#times_scripps <- as.numeric(metadata_scripps$times)

# create distance matrix
r_hellinger_scripps <- r_hellinger_scripps[,-1]
distance_matrix_scripps <- vegdist(r_hellinger_scripps, method = "bray") # changed to vegdist

d2 <- distance_matrix_scripps


x_scripps <- defineTrajectories(d2, entities_scripps, surveys_scripps)

class(x_scripps)

# This object contains two elements:

names(x_scripps)

# Element d contains the input distance matrix, whereas metadata is a data frame including information of observations:

x_scripps$metadata  



trajectoryPCoA(x_scripps, traj.colors = rainbow(n = 14), lwd = 2,
               survey.labels = T)
legend("topright", col= rainbow(n = 14), 
       legend=as.character(c(2011:2020, 2022:2024)), bty="n", lty=1, lwd = 2)

dev.off()

# for(y in 2011:2024){
for(y in c(2011:2020, 2022:2024)){

  
  png(filename = paste0("Figures/ETA/2026-04-18_ETA_", y, ".png"), width = 700, height = 500)
  
  my.ETA <- trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = as.character(y)), 
                 traj.colors = rainbow(n = 14)[y-2010], lwd = 2,
                 survey.labels = T)
  legend("topright", col=rainbow(n = 14)[y-2010], 
         legend=as.character(y), bty="n", lty=1, lwd = 2)
  
  print(my.ETA)
  
  dev.off()
  
  
}




# ######2018 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2018")), 
#                traj.colors = c("#BDB216"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#BDB216"), 
#        legend=c("2018"), bty="n", lty=1, lwd = 2)
# 
# ##### 2019 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2019")), 
#                traj.colors = c("#0072B2"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#0072B2"), 
#        legend=c("2019"), bty="n", lty=1, lwd = 2)
# 
# ##### 2020 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2020")), 
#                traj.colors = c("#CC79A7"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#CC79A7"), 
#        legend=c("2020"), bty="n", lty=1, lwd = 2)
# 
# ###### 2021 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2021")), 
#                traj.colors = c("#56B4E9"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#56B4E9"), 
#        legend=c("2021"), bty="n", lty=1, lwd = 2)
# 
# 
# ###### 2022 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2022")), 
#                traj.colors = c("black"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("black"), 
#        legend=c("2022"), bty="n", lty=1, lwd = 2)
# 
# 
# ##### 2023 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2023")), 
#                traj.colors = c("#D55E00"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#D55E00"), 
#        legend=c("2023"), bty="n", lty=1, lwd = 2)
# 
# ##### 2024 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2024")), 
#                traj.colors = c("#009E73"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#009E73"), 
#        legend=c("2024"), bty="n", lty=1, lwd = 2)


trajectoryLengths(x_scripps)
trajectorySpeeds(x_scripps)
trajectoryAngles(x_scripps)
trajectoryDirectionality(x_scripps)
trajectoryMetrics(x_scripps)
trajectoryInternalVariation(x_scripps)

trajectoryShifts(subsetTrajectories(x_scripps, c("2022","2023")))
trajectoryConvergence(x_scripps, type = "pairwise.asymmetric")
# trajectoryConvergence(x_scripps, type = "multiple")   
#### might have to take the "column" surveys out and instead use the "times" as the survey date 

ScrippsConv.temp <- trajectoryConvergence(x_scripps, type = "pairwise.asymmetric")

corrplot(matrix(as.vector(ScrippsConv.temp$tau)*as.numeric(ScrippsConv.temp$p.value<0.05),13,13))




###############################################################
#################SCRIPPS 2018 to 2024 ETA (only closest days)
###############################################################

# make a data frame from the sample_data
# samples_df_scripps <- data.frame(sample_data(r_hellinger_scripps))

#assigning year and month
# samples_df_scripps$day<-day(samples_df_scripps$Date)  ###assigning day of the year


# STEP 1: Prepare data (use your actual data frame name)
df <- samples_df_scripps %>%
  tibble::rownames_to_column("rownames") %>%
  mutate(
    year = year(Date),
    day = yday(Date)  # use day-of-year for matching
  )

# STEP 2: Define reference year and get all other years
ref_year <- 2011   ##### Isi used 2024 as year with less samples, I am using 2011 as first year in time series
ref_df <- df %>% filter(year == ref_year)
other_years <- setdiff(unique(df$year), ref_year)

# STEP 3: Pairwise matching to closest date within ±5 days (change if needed)
final_matches <- list()

for (yr in other_years) {
  other_df <- df %>% filter(year == yr)
  
  df_ref <- df %>%
    dplyr::rename(ref = rownames, ref_date = Date, ref_day = day)
  
  df_other <- df %>%
    dplyr::rename(other = rownames, other_date = Date, other_day = day)
  
  matched_pairs <- expand.grid(
    ref = ref_df$rownames,
    other = other_df$rownames,
    stringsAsFactors = FALSE
  ) %>%
    left_join(df_ref, by = "ref") %>%
    left_join(df_other, by = "other") %>%
    mutate(date_diff = abs(ref_day - other_day)) %>%
    filter(date_diff <= 5) %>%
    arrange(date_diff)
  
  # Greedy 1-to-1 matching
  used_ref <- character(0)
  used_other <- character(0)
  final_pairs <- list()
  
  for (i in seq_len(nrow(matched_pairs))) {
    row <- matched_pairs[i, ]
    if (!(row$ref %in% used_ref) && !(row$other %in% used_other)) {
      final_pairs[[length(final_pairs) + 1]] <- row
      used_ref <- c(used_ref, row$ref)
      used_other <- c(used_other, row$other)
    }
  }
  
  final_matches[[as.character(yr)]] <- bind_rows(final_pairs)
}

# STEP 4: Keep only reference samples matched in ALL other years
common_refs <- Reduce(intersect, lapply(final_matches, \(x) x$ref))

# STEP 5: Filter all matches to those common refs only
filtered_matches <- lapply(final_matches, \(x) x %>% filter(ref %in% common_refs))

# STEP 6: Gather ref (2024) and matched other years
ref_rows_df <- df %>% filter(rownames %in% common_refs)
other_rows_df <- filtered_matches %>%
  lapply(\(x) df %>% filter(rownames %in% x$other)) %>%
  bind_rows()

# STEP 7: Combine into final matched data frame
matched_df <- bind_rows(ref_rows_df, other_rows_df) %>%
  column_to_rownames("rownames")

# Optional: check balance across years
table(matched_df$year)



matched_df <- matched_df %>%
  rownames_to_column("rownames") %>%
  arrange(year, day) %>%
  column_to_rownames("rownames")


matched_df$point<- "1"


# matched_df %>% count(year)


samples_scripps <- ggplot(matched_df, aes(x = day, y = point)) +
  facet_grid(rows = vars(year), scales = "free", space= "free_x") +
  geom_point(color = "#009E73") +
  theme(axis.text.x = element_text(face="plain", color="black", 
                                   size=5, angle=45, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", 
                                   size=10, angle=0),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"))+
  ylab("sample collection") +
  ggtitle("Scripps Pier time series (closest days ~5)") 
samples_scripps

ggsave(file="samples_scripps_matched.png", samples_scripps, path = "Figures/", width =20 , height = 18, units = "cm")





# Extract sample names (rownames) to keep
samples_to_keep <- rownames(matched_df)

# Subset phyloseq object (replace `r_hellinger_scripps` with your actual object)
ps_matched <- r_hellinger_scripps[which(rownames(r_hellinger_scripps) %in% samples_to_keep),]

# # Prune taxa not present in the subset
# ps_matched <- prune_taxa(taxa_sums(ps_matched) > 0, ps_matched)




metadata_scripps <- matched_df[, c("year", "day")]

# Let us first define the vectors that describe the ecological entity and the survey of each observation:

names(metadata_scripps) <- c("entities", "surveys")

entities_scripps <- as.character(metadata_scripps$entities)
surveys_scripps <- as.integer(metadata_scripps$surveys)
#times_scripps <- as.numeric(metadata_scripps$times)

# create distance matrix
distance_matrix_scripps <- vegdist(ps_matched, method = "bray") # used vegdist here

d2 <- distance_matrix_scripps


x_scripps <- defineTrajectories(d2, entities_scripps, surveys_scripps)

class(x_scripps)

#This object contains two elements:

names(x_scripps)

#Element d contains the input distance matrix, whereas metadata is a data frame including information of observations:

x_scripps$metadata  

dev.off()

for(y in c(2011:2020, 2022:2024)){
  
  
  png(filename = paste0("Figures/ETA/2026-04-18_ETA2_", y, ".png"), width = 700, height = 500)
  
  my.ETA <- trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = as.character(y)), 
                           traj.colors = rainbow(n = 14)[y-2010], lwd = 2,
                           survey.labels = T)
  legend("topright", col=as.character(y), 
         legend=as.character(y), bty="n", lty=1, lwd = 2)
  
  print(my.ETA)
  
  dev.off()
  
  
}




# ######2018 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2018")), 
#                traj.colors = c("#BDB216"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#BDB216"), 
#        legend=c("2018"), bty="n", lty=1, lwd = 2)
# 
# ##### 2019 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2019")), 
#                traj.colors = c("#0072B2"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#0072B2"), 
#        legend=c("2019"), bty="n", lty=1, lwd = 2)
# 
# ##### 2020 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2020")), 
#                traj.colors = c("#CC79A7"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#CC79A7"), 
#        legend=c("2020"), bty="n", lty=1, lwd = 2)
# 
# ###### 2021 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2021")), 
#                traj.colors = c("#56B4E9"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#56B4E9"), 
#        legend=c("2021"), bty="n", lty=1, lwd = 2)
# 
# 
# ###### 2022 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2022")), 
#                traj.colors = c("black"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("black"), 
#        legend=c("2022"), bty="n", lty=1, lwd = 2)
# 
# 
# ##### 2023 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2023")), 
#                traj.colors = c("#D55E00"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#D55E00"), 
#        legend=c("2023"), bty="n", lty=1, lwd = 2)
# 
# ##### 2024 scripps
# trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2024")), 
#                traj.colors = c("#009E73"), lwd = 2,
#                survey.labels = T)
# legend("topright", col=c("#009E73"), 
#        legend=c("2024"), bty="n", lty=1, lwd = 2)



# setwd("~/Documents/Post-doc/03 Bowman Lab/data/20250420_18S/ETA/Figures")


# Open a pdf file
pdf("Figures/ETA/ETA_closest_scripps2011-2024.pdf", width = 25, height = 20) 
# 2. Create a plot
par(mfrow=c(3,3),cex.lab = 3, mar = c(5, 6, 4, 3))   # More space between columns using mar
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2011")), 
               traj.colors = rainbow(n = 14)[1], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2012")), 
               traj.colors = rainbow(n = 14)[2], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2013")), 
               traj.colors = rainbow(n = 14)[3], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2014")), 
               traj.colors = rainbow(n = 14)[4], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2015")), 
               traj.colors = rainbow(n = 14)[5], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2016")), 
               traj.colors = rainbow(n = 14)[6], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2017")), 
               traj.colors = rainbow(n = 14)[7], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2018")), 
               traj.colors = rainbow(n = 14)[8], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2019")), 
               traj.colors = rainbow(n = 14)[9], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2020")), 
               traj.colors = rainbow(n = 14)[10], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2021")), 
               traj.colors = rainbow(n = 14)[11], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2022")), 
               traj.colors = rainbow(n = 14)[12], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2023")), 
               traj.colors = rainbow(n = 14)[13], lwd = 3,
               survey.labels = T)
trajectoryPCoA(subsetTrajectories(x_scripps, site_selection = c("2024")), 
               traj.colors = rainbow(n = 14)[14], lwd = 3,
               survey.labels = T)
# Close the pdf file
dev.off() 


trajectoryLengths(x_scripps)
trajectorySpeeds(x_scripps)
trajectoryAngles(x_scripps)
trajectoryDirectionality(x_scripps)
trajectoryMetrics(x_scripps)
trajectoryInternalVariation(x_scripps)

# trajectoryShifts(subsetTrajectories(x_scripps, c("2022","2023")))

ScrippsConv <- trajectoryConvergence(x_scripps, type = "pairwise.symmetric")
#trajectoryConvergence(x_scripps, type = "multiple")   


library(corrplot)

corrplot(matrix(as.vector(ScrippsConv$tau)*as.numeric(ScrippsConv$p.value<0.05),13,13))

# Open a pdf file
pdf("Figures/ETA/ETA_convergence_scripps2018-2024.pdf", width = 10, height = 10) 
# 2. Create a plot
corrplot(matrix(as.vector(ScrippsConv$tau)*as.numeric(ScrippsConv$p.value<0.05),13,13))
dev.off()



scripps_D_traj_man <- trajectoryDistances(x_scripps, distance.type="DSPD")
print(round(scripps_D_traj_man,3))
scripps_D_traj_man_mat <- as.matrix(scripps_D_traj_man)
scripps_D_traj_man_mat[scripps_D_traj_man_mat == 0] <- NA
corrplot(scripps_D_traj_man_mat, is.corr = F, col = viridis(100))


trajectoryConvergencePlot(x_scripps, type = "pairwise.symmetric")











# setwd("~/Documents/Post-doc/03 Bowman Lab/data/20250420_18S/ETA/Figures")

# Open a pdf file
pdf("Figures/ETA/ETA_closest_scripps2018-2024_all_onebyone.pdf", width = 25, height = 20) 
# 2. Create a plot
par(mfrow=c(3,3),cex.lab = 3, mar = c(5, 6, 4, 3))   # More space between columns using mar
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1], lwd = 3,
               survey.labels = F)
legend("bottomright", col=rainbow(n=14),
       legend=as.character(c(2011:2020, 2022:2024)), bty="n", lty=1, lwd = 3, cex = 2.5)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:2], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:3], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:4], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:5], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:6], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:7], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:8], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:9], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:10], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:11], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:12], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:13], lwd = 3,
               survey.labels = F)
trajectoryPCoA(x_scripps, traj.colors = rainbow(n=14)[1:14], lwd = 3,
               survey.labels = F)
# Close the pdf file
dev.off() 




