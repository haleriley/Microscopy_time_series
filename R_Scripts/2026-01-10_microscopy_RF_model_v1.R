# random forest model of microscropy time series
# RJH
# 2026-01-10

# ---- library ----

library(tidyverse)
library(lubridate)
library(janitor)
library(vegan)
library(plotly)
library(ranger)
library(readxl)
library(patchwork)
library(scico)


# ---- read in data ----

setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/")

weallen.diatoms <- read.csv("AllenData/SIOdiat.csv")
weallen.dinos <- read.csv("AllenData/SIOdino.csv", check.names = F)

ifcb <- read.csv("R_Data/spc_hab_daily_2017.csv")

o2bio <- readRDS("../O2-Ar_time_series/R_Data/2025-08-28_o2bio_est_df.rds")
o2bio.daily <- readRDS("../O2-Ar_time_series/R_Data/2025-10-07_o2bio_est_df_daily.rds")

# hab <- read.csv("https://erddap.sccoos.org/erddap/tabledap/HABs-ScrippsPier.csv?Location_Code%2Ctime%2CTemp%2CAir_Temp%2CSalinity%2CAvg_Chloro%2CAvg_Phaeo%2CPhosphate%2CSilicate%2CNitrite%2CNitrite_Nitrate%2CAmmonium%2CNitrate%2CDA_Volume_Filtered%2CpDA%2CtDA%2CdDA%2CVolume_Settled_for_Counting%2CAkashiwo_sanguinea%2CAlexandrium_spp%2CDinophysis_spp%2CLingulodinium_polyedra%2CProrocentrum_spp%2CPseudo_nitzschia_delicatissima_group%2CPseudo_nitzschia_seriata_group%2CCeratium_spp%2CCochlodinium_spp%2CGymnodinium_spp%2COther_Diatoms%2COther_Dinoflagellates%2CTotal_Phytoplankton&Location_Code=%22SIO%22&time%3E=2008-06-30T15%3A00%3A00Z&time%3C=2025-12-08T20%3A17%3A00Z")
  
sccoos.phytos <- read_excel("R_Data/Microscopy_Carter/SIO_complete_2008chl_2008-2024cells.xlsx")

unique.dinos <- readRDS("R_Data/2026-02-09_unique_dinos.rds")
unique.diatoms <- readRDS("R_Data/2026-02-09_unique_diatomss.rds")

# # ---- clean and reformat ifcb time series data ----
# 
# ifcb <- as.data.frame(ifcb)
# colnames(ifcb)[1] <- "Date"
# ifcb$Date <- parse_date_time(ifcb$Date, orders = "mdY")
# my.dates <- ifcb$Date
# 
# ifcb <- ifcb %>% mutate_all(~ifelse(is.nan(.), NA, .))
# ifcb$Date <- my.dates
# 
# ifcb <- na.omit(ifcb)

# # ---- clean and reformat HAB microscopy time series data ----
# 
# hab.all <- hab[-1,] 
# 
# hab.taxa <- colnames(hab)[19:ncol(hab)]
# 
# hab <- hab[-1,c(2,19:ncol(hab))]
# colnames(hab)[1] <- "Date"
# 
# hab$Date <- parse_date_time(hab$Date, orders = "Ymd HMS")
# hab$Date <- parse_date_time(paste(year(hab$Date), month(hab$Date), day(hab$Date), sep = "-"), orders = "Ymd")
# my.dates <- hab$Date
# 
# hab <- hab %>% mutate(across(where(is.character), ~ na_if(.x, "NaN")))
# 
# hab$Date <- my.dates
# hab <- na.omit(hab)
# 

# ---- clean and reformat phyto count data ----

sccoos.phytos <- sccoos.phytos[,c(1,21:153)]

sccoos.phytos <- clean_names(sccoos.phytos)

sccoos.phytos[] <- lapply(sccoos.phytos, as.numeric)

sccoos.phytos <- sccoos.phytos[which(sccoos.phytos$total_phytoplankton > 0),]

test <- t(sccoos.phytos[146:156,])
testtest <- as.data.frame(rowSums(test))
## nitzschia_spp are not routinely collected throughout, so removing here
sccoos.phytos <- sccoos.phytos[,-which(colnames(sccoos.phytos) == "nitzschia_spp")]

## total zooplankton and "other" counts missing for 2015, not really important for this model so removing that column too
zoop.other.col.index <- c(118:132)
sccoos.phytos <- sccoos.phytos[,-zoop.other.col.index]

sccoos.phytos <- na.omit(sccoos.phytos)


sccoos.phytos$sample_id <- parse_date_time(sccoos.phytos$sample_id, orders = "ymd")
colnames(sccoos.phytos)[1] <- "Date"

overlapping.taxa <- colnames(sccoos.phytos)[which(!(colnames(sccoos.phytos) %in% c("Date", "total_diatoms", "total_dinoflagellates", "total_phytoplankton")))]

# # ---- clean W.E. Allen dino and diatom data frames ----
# 
# weallen.diatoms <- clean_names(weallen.diatoms)
# weallen.dinos <- clean_names(weallen.dinos) # ignore warnings...I think....
# 
# weallen.diatoms$Date <- parse_date_time(paste(weallen.diatoms$year,weallen.diatoms$month, weallen.diatoms$day), orders = "Ymd")
# weallen.diatoms <- weallen.diatoms[,-c(1:3)]
# 
# weallen.dinos$Date <- parse_date_time(paste(weallen.dinos$year,weallen.dinos$month, weallen.dinos$day), orders = "Ymd")
# weallen.dinos <- weallen.dinos[,-c(1:3)]
# 
# # ---- split SCCOOS phytos into diatoms and dinoflagellates to match how relative abundance is calculated in W.E. Allen time series ----
# 
# sccoos.diatoms <- sccoos.phytos[,c(1,2:62)]
# sccoos.dinos <- sccoos.phytos[,c(1,66:114)]
# 
# 
# # ---- change absolute counts to relative abundance ----
# 
# ## changing absolute counts to relative abundance to more easily match 
# ## W.E. Allen time series, since methods have changed
# 
# ## need to do this for diatoms and dinoflagellates separately 
# 
# sccoos.dinos.dates <- sccoos.dinos$Date
# sccoos.dinos <- sccoos.dinos[,-1]
# sccoos.dinos <- sccoos.dinos/rowSums(sccoos.dinos)
# sccoos.dinos$Date <- sccoos.dinos.dates
# 
# sccoos.diatoms.dates <- sccoos.diatoms$Date
# sccoos.diatoms <- sccoos.diatoms[,-1]
# sccoos.diatoms <- sccoos.diatoms/rowSums(sccoos.diatoms)
# sccoos.diatoms$Date <- sccoos.diatoms.dates
# 
# weallen.dinos[is.na(weallen.dinos)] <- 0
# weallen.dinos.dates <- weallen.dinos$Date
# weallen.dinos <- weallen.dinos[,-which(colnames(weallen.dinos) == "Date")]
# weallen.dinos <- weallen.dinos/rowSums(weallen.dinos)
# weallen.dinos$Date <- weallen.dinos.dates
# 
# weallen.diatoms[is.na(weallen.diatoms)] <- 0
# weallen.diatoms.dates <- weallen.diatoms$Date
# weallen.diatoms <- weallen.diatoms[,-which(colnames(weallen.diatoms) == "Date")]
# weallen.diatoms <- weallen.diatoms/rowSums(weallen.diatoms)
# weallen.diatoms$Date <- weallen.diatoms.dates
# 
# # ---- combine dino and diatom dataframes ----
# 
# weallen.combo <- merge(weallen.diatoms, weallen.dinos, by = "Date")
# 
# sccoos.combo <- merge(sccoos.diatoms, sccoos.dinos, by = "Date")
# 

# # ---- manually correct taxonomy names to match using combine_phyto_names.xlsx ----
# 
# # overlappusing ing.taxa.index <- which(colnames(weallen.combo)[-1] %in% colnames(phytos)[-1])
# 
# # my.df1 <- as.data.frame(colnames(weallen.combo)[-1])
# # my.df2 <- as.data.frame(colnames(phytos)[-1])
# # my.df <- as.data.frame(colnames(weallen.combo)[overlapping.taxa.index])
# 
# combine.taxonomies.df <- read_excel("combine_phyto_names.xlsx")
# 
# new.tax.weallen <- combine.taxonomies.df[,c(1,4)]
# index <- which(is.na(new.tax.weallen$New_name_WEAllen) == F)
# new.tax.weallen$WEAllen[index] <- new.tax.weallen$New_name_WEAllen[index]
# 
# new.tax.sccoos <- combine.taxonomies.df[,c(5,6)]
# index <- which(is.na(new.tax.sccoos$New_name_SCCOOS) == F)
# new.tax.sccoos$SCCOOS[index] <- new.tax.sccoos$New_name_SCCOOS[index]
# 
# index <- which(is.na(new.tax.sccoos$SCCOOS) == F)
# new.tax.sccoos <- new.tax.sccoos[index,]
# colnames(phytos)[-1] <- new.tax.sccoos$SCCOOS
# 
# colnames(weallen.combo)[-1] <- new.tax.weallen$WEAllen
# 
# 
# 
# weallen.combo.long <- weallen.combo %>% pivot_longer(cols = colnames(weallen.combo)[-1], names_to = "Taxa", values_to = "Count")
# weallen.combo.long <- weallen.combo.long %>% group_by(Date,Taxa) %>% summarize(Count = sum(Count, na.rm = T))
# weallen.combo.clean <- weallen.combo.long %>% pivot_wider(id_cols = Date, names_from = Taxa, values_from = Count)
# 
# 
# phytos.long <- phytos %>% pivot_longer(cols = colnames(phytos)[-1], names_to = "Taxa", values_to = "Count")
# phytos.long <- phytos.long %>% group_by(Date,Taxa) %>% summarize(Count = sum(Count, na.rm = T))
# phytos.clean <- phytos.long %>% pivot_wider(id_cols = Date, names_from = Taxa, values_from = Count)
# 

# ---- combine phyto names, excluding spp. groups ----

# overlapping.taxa.index <- which(colnames(weallen.combo)[-1] %in% colnames(sccoos.combo)[-1])
# overlapping.taxa <- colnames(weallen.combo)[-1][which(colnames(weallen.combo)[-1] %in% colnames(sccoos.combo)[-1])]
# 
# weallen.combo <- weallen.combo[,which(colnames(weallen.combo) %in% c("Date", overlapping.taxa))]
sccoos.combo <- sccoos.phytos[,which(colnames(sccoos.phytos) %in% c("Date", overlapping.taxa))]


# weallen.combo.long <- weallen.combo %>% pivot_longer(cols = colnames(weallen.combo)[-1], names_to = "Taxa", values_to = "Count")
# weallen.combo.long <- weallen.combo.long %>% group_by(Date,Taxa) %>% summarize(Count = sum(Count, na.rm = T))
# weallen.combo.clean <- weallen.combo.long %>% pivot_wider(id_cols = Date, names_from = Taxa, values_from = Count)


sccoos.combo.long <- sccoos.combo %>% pivot_longer(cols = colnames(sccoos.combo)[-1], names_to = "Taxa", values_to = "Count")
sccoos.combo.long <- sccoos.combo.long %>% group_by(Date,Taxa) %>% dplyr::summarize(Count = sum(Count, na.rm = T))
sccoos.combo.clean <- sccoos.combo.long %>% pivot_wider(id_cols = Date, names_from = Taxa, values_from = Count)



# ---- combine W.E. Allen and SCCOOS phyto time series!!!! ----

# overlapping.taxa <- colnames(weallen.combo.clean)[colnames(weallen.combo.clean) %in% colnames(sccoos.combo.clean)]
# overlapping.taxa <- colnames(sccoos.combo.clean)[colnames(sccoos.combo.clean) %in% colnames(weallen.combo.clean)]
# 
# sccoos.combo.clean <- sccoos.combo.clean[,overlapping.taxa]
# weallen.combo.clean <- weallen.combo.clean[,overlapping.taxa]
# 
# phytos.full <- rbind(weallen.combo.clean, sccoos.combo.clean)
# phytos.full <- phytos.full[,c(ncol(phytos.full),2:ncol(phytos.full)-1)]
phytos.full <- sccoos.combo.clean

# # ---- optional: hellinger transformation of phyto rel abunds ----
# 
# phytos.full <- na.omit(phytos.full)
# phytos.dates <- phytos.full$Date
# phytos.full <- phytos.full[,-which(colnames(phytos.full) == "Date")]
# phytos.full <- decostand(phytos.full, method = "hellinger")
# phytos.full$Date <- phytos.dates

# ---- smooth/average o2bio data ----

o2bio.smooth <- o2bio[1,]

combo <- merge(o2bio.daily, phytos.full, by = "Date")


for(d in 1:nrow(combo)){
  
  my.date <- combo$Date[d]
  
  my.date1 <- my.date - 60*60*24
  my.date2 <- my.date + 60*60*47
  
  my.df <- o2bio[which(o2bio$Date.Time >= my.date1 & o2bio$Date.Time <= my.date2),]
  
  my.df <- my.df %>% summarize_all(mean)
  
  my.df$Date.Time <- my.date
  
  o2bio.smooth <- rbind(o2bio.smooth, my.df)
  
}

o2bio.smooth <- o2bio.smooth[-1,]

colnames(o2bio.smooth)[which(colnames(o2bio.smooth) == "Date.Time")] <- "Date"



combo <- merge(o2bio.smooth[which(colnames(o2bio.smooth) %in% c("Date", "O2bio.estimated"))], phytos.full, by = "Date")

predictors <- colnames(phytos.full)[-1]
response <- "O2bio.estimated"

saveRDS(predictors, file = "2026-03-27_phyto_rf_predictors.rds")

# ---- split training and testing data ----


0.1*difftime(min(combo$Date), max(combo$Date), units = "days")
0.1*nrow(combo)
test.index <- sample(c(1:nrow(combo)), size = 30, replace = F)
mean(difftime(combo$Date[test.index], combo$Date[test.index + 21], units = "days"), na.rm = T)

date1 <- parse_date_time("2019-08-01", orders = "Ymd")
date2 <- parse_date_time("2019-11-01", orders = "Ymd")

test.index <- which(combo$Date >= date1 & combo$Date <= date2)
train <- combo[-test.index,]
test <- combo[test.index,]

nrow(test)/nrow(combo)

ggplot(combo) +
  geom_rect(aes(xmin = date1, xmax = date2, ymin = -Inf, ymax = Inf), color = "blue") +
  geom_line(aes(x = Date, y = O2bio.estimated)) +
  theme_bw()


# ---- o2bio prediction model ----

set.seed(1234)
m1 <- ranger(as.formula(paste(response, '.', sep = '~')),
             data = train[,c(response, predictors)])

## two different ways of calculating RMSE
sqrt(mean((m1$predictions - train$O2bio.estimated)^2))
m1.RMSE <- sqrt(m1$prediction.error)

plot(m1$predictions ~ train[,response])

## assess model performance by testing model on witheld data (test data)
set.seed(1234)
o2bio.predict <- predict(m1, test) 

plot(o2bio.predict$predictions ~ test[,response],
     ylab = 'Observed',
     xlab = 'Predicted')

m1.lm <- lm(o2bio.predict$predictions ~ test[,response])
abline(0, 1, lty = 2)
abline(m1.lm)
summary(m1.lm)

# ---- parameter optimization ----

## define the parameter space
## these are all of the little settings in the model that we will try and adjust to find the best fit for the data

hyper.grid <- expand.grid(
  n.edges = seq(100, 3000, 100), # aka n.trees
  mtry       = seq(1, ncol(ifcb), by = 1), # 
  node_size  = seq(3, 9, by = 2), 
  sample_size = c(.55, .632, .70, .80), # internal
  OOB_RMSE   = 0 # internal
)

set.seed(1234)
for(i in 1:nrow(hyper.grid)){ ## AKA for every combination of parameter settings
  
  # predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)] # cannot use more n.edges than boruta predictors
  
  try({ ## try clause necessary because some parameter combinations are incompatible
    
    model <- ranger(
      formula = as.formula(paste(response, '.', sep = '~')),
      data = train[,c(response, predictors)], 
      num.trees       = 500,
      mtry            = hyper.grid$mtry[i],
      min.node.size   = hyper.grid$node_size[i],
      sample.fraction = hyper.grid$sample_size[i],
      seed            = 123
    )
    
    ## add OOB error to grid
    hyper.grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
    
    ## From the internet: 
    ## OOB (out-of-bag) score is a performance metric for a machine learning model, 
    ## specifically for ensemble models such as random forests. 
    ## It is calculated using the samples that are not used in the training of the model, 
    ## which is called out-of-bag samples.
    ## The OOB_score is computed as the number of correctly predicted rows from the out-of-bag sample. 
    ## OOB Error is the number of wrongly classifying the OOB Sample.
    
  }, silent = F)
  
  print(paste(i, 'out of', nrow(hyper.grid), hyper.grid$OOB_RMSE[i]))
  
}

hyper.grid$OOB_RMSE[hyper.grid$OOB_RMSE == 0] <- NA
hyper.grid <- na.omit(hyper.grid)

hist(hyper.grid$OOB_RMSE, breaks = 100)

## define selected optimal parameters for the model
selected.params <- hyper.grid[which.min(hyper.grid$OOB_RMSE),]

saveRDS(selected.params, "2026-04-26_phyto_rf_selected_params.rds")

# ---- apply hypertuned parameters to model ----

# predictors <- boruta.index[order(colSums(asv.train)[colnames(asv.train) %in% boruta.index], decreasing = T)]

## create second model using optimal selected paramters
set.seed(1234)
m2 <- ranger(
  formula = as.formula(paste(response, '.', sep = '~')),
  data = train[,c(response, predictors)],
  num.trees       = 500,
  mtry            = selected.params$mtry,
  min.node.size   = selected.params$node_size,
  sample.fraction = selected.params$sample_size,
  seed            = 123,
  importance = 'permutation',
  oob.error = T
)



# ---- create final model ----

set.seed(1234)
## create second model using optimal selected parameters
m3 <- ranger(
  formula = as.formula(paste(response, '.', sep = '~')),
  data = combo[,c(response, predictors)],
  num.trees       = 500,
  mtry            = selected.params$mtry,
  min.node.size   = selected.params$node_size,
  sample.fraction = selected.params$sample_size,
  seed            = 123,
  importance = 'permutation',
  oob.error = T
)


set.seed(1234)
## compare m1 and m2 with and without parameter optimization
o2bio.predict <- predict(m1, test)
o2bio.predict <- predict(m2, test)
o2bio.predict <- predict(m3, test)

m3.RMSE <- mean(sqrt((test$O2bio.estimated - o2bio.predict$predictions)^2))
m3.rel.error <- median(abs(test$O2bio.estimated - o2bio.predict$predictions)/abs(test$O2bio.estimated))
summary(lm(o2bio.predict$predictions~test$O2bio.estimated))

o2bio.predict <- predict(m3, combo)

plot(o2bio.predict$predictions ~ combo[,response],
     ylab = 'Predicted',
     xlab = 'Observed')
abline(0, 1, lty = 2)
abline(lm(o2bio.predict$predictions ~ combo[,response]))
summary(lm(o2bio.predict$predictions ~ combo[,response]))

set.seed(1234)
o2bio.predict <- predict(m3, combo)

saveRDS(o2bio.predict, "202-04-27_o2bio_predict_m3.rds")


# mean(sqrt((train$o2bio - o2bio.predict$predictions)^2))
# median(abs(train$o2bio - o2bio.predict$predictions)/abs(train$o2bio))

combo$O2bio.predicted <- o2bio.predict$predictions

m3.RMSE.full <- mean(sqrt((combo$O2bio.estimated - combo$O2bio.predicted)^2))
m3.rel.error.full <- median(abs(combo$O2bio.estimated - combo$O2bio.predicted)/abs(combo$O2bio.estimated))



# ---- plot time series ----

ggplot(data = combo) +
  geom_line(aes(x = Date, y = O2bio.estimated), color = "violetred1", linewidth = 1) +
  geom_line(aes(x = Date, y = O2bio.predicted), color = "violetred4", linewidth = 1) +
  theme_bw() +
  labs(x = "Date", y = "[O2]bio (uM)") 

ggplot(data = combo) +
  geom_abline(aes(slope = 1, intercept = 0), color = "violetred1", linewidth = 1) +
  geom_point(aes(x = O2bio.estimated, y = O2bio.predicted), color = "violetred4") +
  theme_bw() +
  labs(x = "[O2]bio", y = "Predicted [O2]bio")

summary(lm(combo$O2bio.predicted~combo$O2bio.estimated))
my.slope <- summary(lm(combo$O2bio.predicted~combo$O2bio.estimated))$coefficients[2,1]
my.intercept <- summary(lm(combo$O2bio.predicted~combo$O2bio.estimated))$coefficients[1,1]


combo$O2bio.predicted.altered <- ((combo$O2bio.predicted-my.intercept)/my.slope)

ggplot(data = combo) +
  geom_line(aes(x = Date, y = O2bio.estimated), color = "violetred1", linewidth = 1) +
  geom_line(aes(x = Date, y = O2bio.predicted.altered), color = "violetred4", linewidth = 1) +
  theme_bw() +
  labs(x = "Date", y = "[O2]bio (uM)") 

ggplot(data = combo) +
  geom_abline(aes(slope = 1, intercept = 0), color = "violetred1", linewidth = 1) +
  geom_point(aes(x = O2bio.estimated, y = O2bio.predicted.altered), color = "violetred4") +
  theme_bw() +
  labs(x = "[O2]bio", y = "Predicted [O2]bio")

summary(lm(combo$O2bio.predicted.altered~combo$O2bio.estimated))


m3.RMSE.full.altered <- mean(sqrt((combo$O2bio.estimated - combo$O2bio.predicted.altered)^2))
m3.rel.error.full.altered <- median(abs(combo$O2bio.estimated - combo$O2bio.predicted.altered)/abs(combo$O2bio.estimated))


# ---- top predictors ----

model.var.imp <- as.data.frame(m3$variable.importance)
colnames(model.var.imp)[1] <- "Variable.Importance"
model.var.imp$Predictors <- rownames(model.var.imp)

model.var.imp <- model.var.imp[order(model.var.imp$Variable.Importance, decreasing = T),]


unique.diatoms <- unique(colnames(sccoos.phytos[,2:64]))
unique.dinos <- unique(colnames(sccoos.phytos[66:116]))
unique.diatoms <- na.omit(unique.diatoms)
unique.dinos <- na.omit(unique.dinos)

saveRDS(unique.dinos, file = "2026-04-23_unique_dinos.rds")
saveRDS(unique.diatoms, file = "2026-04-23_unique_diatomss.rds")


model.var.imp$Group <- NA
model.var.imp$Group[which(model.var.imp$Predictors %in% unique.diatoms)] <- "Diatom"
model.var.imp$Group[which(model.var.imp$Predictors %in% unique.dinos)] <- "Dinoflagellate"

model.var.imp$Predictors <- factor(model.var.imp$Predictors, levels = model.var.imp$Predictors[order(model.var.imp$Variable.Importance, decreasing = F)])

# ggplot(data = model.var.imp) +
#   geom_bar(aes(x = Predictors, y = Variable.Importance, fill = Variable.Importance), stat = "identity") +
#   theme_bw() +
#   scale_fill_viridis_c() +
#   coord_flip() +
#   theme(legend.position = "none") +
#   theme(axis.text.y = element_text(color = rev(model.var.imp$Group))) # not sure why this needs to be reversed but it does
#   
  
model.var.imp <- model.var.imp[which(model.var.imp$Variable.Importance > 0),]

ggplot(data = model.var.imp) +
  geom_bar(aes(x = Predictors, y = Variable.Importance, fill = Group), stat = "identity") +
  theme_bw() +
  # scale_fill_viridis_c() +
  scale_fill_manual(values = c("deepskyblue2", "deepskyblue4")) +
  labs(y = "Variable Model Influence") +
  coord_flip() +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  labs(tag = "A") + theme(plot.tag = element_text(size = 18, face = "bold")) 

saveRDS(model.var.imp, file = "2026-04-26_microscopy_model_var_imp.rds")

saveRDS(predictors, file = "2026-04-2_microscopy_model_predictors.rds")

# ---- apply model to full phyto time series ----

# o2bio.predict <- predict(m3, combo)

combo.full <- merge(o2bio.smooth[which(colnames(o2bio.smooth) %in% c("Date", "O2bio.estimated"))], phytos.full, by = "Date", all = T)

full.phyto.o2bio.predict <- predict(m3, combo.full)

combo.full$O2bio.predicted <- full.phyto.o2bio.predict$predictions
combo.full$O2bio.predicted.altered <- (full.phyto.o2bio.predict$predictions-my.intercept)/my.slope

combo.full$Dataset <- "W.E. Allen 2"
combo.full$Dataset[which(combo.full$Date <= parse_date_time("1929", orders = "Y"))] <- "W.E. Allen 1"
combo.full$Dataset[which(combo.full$Date >= parse_date_time("2000", orders = "Y"))] <- "SCCOOS"
  

# saveRDS(combo.full, file = "2026-02-09_microscopy_o2bio_combo_full.rds")
# 
# combo.full <- readRDS("2026-02-09_microscopy_o2bio_combo_full.rds")


saveRDS(combo.full, file = "2026-04-23_microscopy_o2bio_combo_full.rds")



setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/R_Data/")

combo.full <- readRDS("2026-04-23_microscopy_o2bio_combo_full.rds")


my.min <- min(na.omit(c(combo.full$O2bio.estimated, combo.full$O2bio.predicted.altered)))
my.max <- max(na.omit(c(combo.full$O2bio.estimated, combo.full$O2bio.predicted.altered)))

n.years.sccoos <- length(unique(year(combo.full$Date[which(combo.full$Dataset == "SCCOOS")])))
n.years.weallen <- length(unique(year(combo.full$Date[which(combo.full$Dataset == "W.E. Allen 1" | combo.full$Dataset == "W.E. Allen 2")])))
n.years.total <- sum(n.years.sccoos, n.years.weallen)

# a <- ggplot() +
#   geom_hline(yintercept = 0, alpha = 0.5) +
#   # geom_line(aes(x = Date, y = O2bio.estimated), color = "violetred1", linewidth = 1) +
#   geom_line(data = combo.full[which(combo.full$Dataset == "W.E. Allen 1"),], aes(x = Date, y = O2bio.predicted.altered), color = "violetred4", linewidth = 1, alpha = 0.7) +
#   geom_line(data = combo.full[which(combo.full$Dataset == "W.E. Allen 2"),], aes(x = Date, y = O2bio.predicted.altered), color = "violetred4", linewidth = 1, alpha = 0.7) +
#   # geom_hline(aes(yintercept = mean(combo.full$O2bio.estimated)), color = violetred4) +
#   geom_hline(aes(yintercept = mean(combo.full$O2bio.predicted.altered[which(combo.full$Dataset == "W.E. Allen 1" | combo.full$Dataset == "W.E. Allen 2")])), color = "violetred4") +
#   theme_bw() +
#   labs(x = "Date", y = "[O2]bio (uM)") +
#   scale_y_continuous(breaks = seq(from = round(my.min, -1), round(my.max, -1), by = 20), limits = c(my.min, my.max)) +
#   # ylim(c(my.min, my.max)) +
#   scale_x_datetime(date_breaks = "3 years", date_labels = "%Y", date_minor_breaks = "1 year")

b <- ggplot(data = combo.full[which(combo.full$Dataset == "SCCOOS"),]) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_line(aes(x = Date, y = O2bio.predicted.altered, color = "Predicted [O2]bio"), linewidth = 1, alpha = 0.7) +
  geom_line(aes(x = Date, y = O2bio.estimated, color = "[O2]bio"), linewidth = 1, alpha = 0.7) +
  geom_hline(aes(yintercept = mean(combo.full$O2bio.predicted.altered[which(combo.full$Dataset == "SCCOOS")])), color = "violetred4") +
  geom_hline(aes(yintercept = mean(combo.full$O2bio.estimated[which(combo.full$Dataset == "SCCOOS")])), color =  "#0A9396") +
  scale_color_manual(name = "", values = c("Predicted [O2]bio" = "violetred4", "[O2]bio" = "#0A9396"), labels = c("[O2]bio", "Predicted [O2]bio")) +
  theme_bw() +
  labs(x = "Date", y = "[O2]bio") +
  scale_y_continuous(breaks = seq(from = round(my.min, -1), round(my.max, -1), by = 20), limits = c(my.min, my.max)) +
  scale_x_datetime(date_breaks = "2 years", date_labels = "%Y", date_minor_breaks = "1 year") +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  labs(tag = "A") + theme(plot.tag = element_text(size = 18, face = "bold")) 

b

# a+b +  plot_layout(widths = c(n.years.weallen/n.years.total, n.years.sccoos/n.years.total))


ggplot(data = combo.full) +
  geom_abline(aes(slope = 1, intercept = 0), color = "#0A9396", linewidth = 1) +
  geom_point(aes(x = O2bio.estimated, y = O2bio.predicted.altered), color = "violetred4", alpha = 0.7) +
  theme_bw() +
  labs(x = "[O2]bio", y = "Predicted [O2]bio") +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  labs(tag = "B") + theme(plot.tag = element_text(size = 18, face = "bold"))

summary(lm(combo.full$O2bio.predicted.altered~combo.full$O2bio.estimated))


m3.RMSE.full.altered <- mean(sqrt((combo.full$O2bio.estimated - combo.full$O2bio.predicted.altered)^2), na.rm = T)
m3.rel.error.full.altered <- median(abs(combo.full$O2bio.estimated - combo.full$O2bio.predicted.altered)/abs(combo.full$O2bio.estimated), na.rm = T)
  




ggplot(data = combo.full[which(combo.full$Dataset == "SCCOOS"),]) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_line(aes(x = Date, y = O2bio.predicted, color = "Predicted [O2]bio"), linewidth = 1, alpha = 0.7) +
  geom_line(aes(x = Date, y = O2bio.estimated, color = "[O2]bio"), linewidth = 1, alpha = 0.7) +
  geom_hline(aes(yintercept = mean(combo.full$O2bio.predicted[which(combo.full$Dataset == "SCCOOS")])), color = "violetred4") +
  geom_hline(aes(yintercept = mean(combo.full$O2bio.estimated[which(combo.full$Dataset == "SCCOOS")])), color =  "#0A9396") +
  scale_color_manual(name = "", values = c("Predicted [O2]bio" = "violetred4", "[O2]bio" = "#0A9396"), labels = c("[O2]bio", "Predicted [O2]bio")) +
  theme_bw() +
  labs(x = "Date", y = "[O2]bio") +
  scale_y_continuous(breaks = seq(from = round(my.min, -1), round(my.max, -1), by = 20), limits = c(my.min, my.max)) +
  scale_x_datetime(date_breaks = "2 years", date_labels = "%Y", date_minor_breaks = "1 year") +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  labs(tag = "A") + theme(plot.tag = element_text(size = 18, face = "bold")) 


# a+b +  plot_layout(widths = c(n.years.weallen/n.years.total, n.years.sccoos/n.years.total))


ggplot(data = combo.full) +
  geom_abline(aes(slope = 1, intercept = 0), color = "#0A9396", linewidth = 1) +
  geom_point(aes(x = O2bio.estimated, y = O2bio.predicted), color = "violetred4", alpha = 0.7) +
  theme_bw() +
  labs(x = "[O2]bio", y = "Predicted [O2]bio") +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  labs(tag = "B") + theme(plot.tag = element_text(size = 18, face = "bold"))

summary(lm(combo.full$O2bio.predicted.altered~combo.full$O2bio.estimated))





# t.test(x = combo.full$O2bio.predicted.altered[which(combo.full$Dataset == "W.E. Allen 1" | combo.full$Dataset == "W.E. Allen 2")],
#        y = combo.full$O2bio.predicted.altered[which(combo.full$Dataset == "SCCOOS")] )


# ---- RF temporal "chunk" cross validation, following Ch. 1 & 2 and Hale et al. (2026) ----

setwd("C://Users/haler/Documents/PhD-Bowman/Microscopy_time_series/R_Data/")

cross.val.df <- na.omit(combo.full)

all.valid.date.times <- c(cross.val.df$Date[1])

0.1*difftime(min(cross.val.df$Date), max(cross.val.df$Date), units = "days")

which(as.character(cross.val.df$Date) == "2021-03-08")

for(d in 1:nrow(cross.val.df)){
  
  date1 <- cross.val.df$Date[d] - 105*60*60*24
  date2 <- cross.val.df$Date[d] + 105*60*60*24  # 90/10 split, ~7 months
  
  my.df <- cross.val.df$Date[which(cross.val.df$Date >= date1 & cross.val.df$Date <= date2)]
  
  if(length(my.df) >= 0.8*30){ # number of expected samples in 7 months *0.8 in case a little missing data
    
    all.valid.date.times <- c(all.valid.date.times, cross.val.df$Date[d]) 
    
  }
  
}

all.valid.date.times <- all.valid.date.times[-1]
all.valid.days <- all.valid.date.times[which(hour(all.valid.date.times) == 0)]
# all.valid.dates <- parse_date_time(paste(year(all.valid.date.times), month(all.valid.date.times), day(all.valid.date.times), sep = "-"), orders = "Ymd")
# all.valid.days <- unique(all.valid.dates)

my.cross.val.values <- as.data.frame(matrix(data = NA, nrow = length(all.valid.days), ncol = 19))
colnames(my.cross.val.values) <- c("Date", "n.days", "perc.data", "RMSE_model_hyper", "R2_model_hyper", "p_model_hyper", "med_perc_error_hyper", "RMSE_model_hyper_full", "R2_model_hyper_full", "p_model_hyper_full", "med_perc_error_hyper_full", "mean.O2bio", "mean.O2bioest", "med.O2bio", "med.O2bioest")

my.cross.val.values$Date <- all.valid.days

selected.params <- readRDS("2026-04-26_phyto_rf_selected_params.rds")

which(as.character(my.cross.val.values$Date) == "2021-03-08")

d=1
for(d in 1:length(all.valid.days)){
  
  date1 <- all.valid.days[d] - 105*60*60*24
  date2 <- all.valid.days[d] + 105*60*60*24  # 90/10 split, ~7 months
  
  index <- which(cross.val.df$Date >= date1 & cross.val.df$Date <= date2)
  
  cross.val.df.nd <- cross.val.df[,which(colnames(cross.val.df) != "Date")]
  
  cross.val.test <- cross.val.df.nd[index,]
  cross.val.train <- cross.val.df.nd[-index,] 
  
  my.cross.val.values$perc.data[d] <- nrow(cross.val.test)/nrow(cross.val.train)
  
  # my.cross.val.values$RMSE_no_model[d] <- sqrt(mean((cross.val.test$O2bio.predicted - cross.val.test$O2bio.estimated)^2))
  # my.cross.val.values$R2_no_model[d] <- summary(lm(formula = cross.val.test$O2bio.predicted~cross.val.test$O2bio.estimated))$r.squared
  # my.cross.val.values$p_no_model[d] <- summary(lm(formula = cross.val.test$O2bio.predicted~cross.val.test$O2bio.estimated))$coefficients[2,4]
  # my.cross.val.values$med_perc_error_no_model[d] <- median((cross.val.test$O2bio.predicted - cross.val.test$O2bio.estimated)/cross.val.test$O2bio.estimated)
  
  temp.rf <- ranger(
    formula = as.formula(paste(response, '.', sep = '~')),
    data = cross.val.train[,c(response, predictors)],
    num.trees       = 500,
    mtry            = selected.params$mtry,
    min.node.size   = selected.params$node_size,
    sample.fraction = selected.params$sample_size,
    seed            = 123,
    importance = 'permutation',
    oob.error = T
  )
  
  
  temp.prediction <- predict(temp.rf, cross.val.test)
  temp.prediction <- temp.prediction$predictions
  
  my.cross.val.values$RMSE_model_hyper[d] <- sqrt(mean((temp.prediction-cross.val.test$O2bio.estimated)^2))
  my.cross.val.values$R2_model_hyper[d] <- summary(lm(formula = temp.prediction~cross.val.test$O2bio.estimated))$r.squared
  my.cross.val.values$p_model_hyper[d] <- summary(lm(formula = temp.prediction~cross.val.test$O2bio.estimated))$coefficients[2,4]
  my.cross.val.values$med_perc_error_hyper[d] <- median((cross.val.test$O2bio.estimated - temp.prediction)/cross.val.test$O2bio.estimated)
  my.cross.val.values$mean.O2bio[d] <- mean(cross.val.test$O2bio.estimated)
  my.cross.val.values$mean.O2biopred[d] <- mean(temp.prediction)
  my.cross.val.values$med.O2bio[d] <- median(cross.val.test$O2bio.estimated)
  my.cross.val.values$med.O2biopred[d] <- median(temp.prediction)
  
  full.prediction <- predict(temp.rf, cross.val.df)
  full.prediction <- full.prediction$predictions
  
  my.cross.val.values$RMSE_model_hyper_full[d] <- sqrt(mean((full.prediction-cross.val.df$O2bio.estimated)^2))
  my.cross.val.values$R2_model_hyper_full[d] <- summary(lm(formula = full.prediction~cross.val.df$O2bio.estimated))$r.squared
  my.cross.val.values$p_model_hyper_full[d] <- summary(lm(formula = full.prediction~cross.val.df$O2bio.estimated))$coefficients[2,4]
  my.cross.val.values$med_perc_error_hyper_full[d] <- median((cross.val.df$O2bio.estimated - full.prediction)/cross.val.df$O2bio.estimated)
  
  
  print(paste("Cross-Validation", d, "of", length(all.valid.days), "complete"))
  
}


saveRDS(my.cross.val.values, file = "2026-04-27_phyto_rf_cross_validation_results.rds")
my.cross.val.values <- readRDS("2026-04-27_phyto_rf_cross_validation_results.rds")

range(my.cross.val.values$perc.data)

# hist(my.cross.val.values$RMSE_no_model)
hist(my.cross.val.values$RMSE_model_hyper)
hist(my.cross.val.values$RMSE_model_hyper_full)

# hist(my.cross.val.values$med_perc_error_no_model)
hist(my.cross.val.values$med_perc_error_hyper)
hist(my.cross.val.values$med_perc_error_hyper_full)

# hist(my.cross.val.values$R2_no_model)
hist(my.cross.val.values$R2_model_hyper)
hist(my.cross.val.values$R2_model_hyper_full)

colnames(my.cross.val.values)[c(7,11)] <- c("2-week Validation", "Full Time Series")
cross.val.long <- my.cross.val.values[,c(1,7)] %>% pivot_longer(cols = colnames(my.cross.val.values)[c(7)], names_to = "model", values_to = "my.perc.error")
cross.val.long$Year <- year(cross.val.long$Date)
cross.val.long$Date <- parse_date_time(paste(month(cross.val.long$Date), day(cross.val.long$Date), sep = "-"), orders = "md")

# ggplot(data = cross.val.long) +
#   geom_point(aes(x = Date, y = abs(my.perc.error)*100, color = factor(model, levels = c("2-week Validation", "Full Time Series"))), size = 2) +
#   facet_wrap(.~Year, ncol = 1) +
#   scale_x_datetime(date_breaks = "1 month", date_labels = "%b") +
#   theme_bw() +
#   labs(x = "Date", y = "Percent error (%)", color = "Model") +
#   scale_color_manual(name = "Model", 
#                      values = c("2-week Validation" = brewer.pal(8,"Dark2")[6], "Full Time Series" = brewer.pal(8,"Dark2")[3]),
#                      labels = c("2-week validation", "Full time series"),
#                      guide = "legend") +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold"))
# 
# ggplot(data = cross.val.long) +
#   geom_line(data = combined.df.overlap, aes(x = Date.Time, y = aou), lwd = 1, alpha = 0, color = "white") +
#   geom_bar(aes(x = Date.Time, y = abs(my.perc.error)*100, fill = factor(model, levels = c("2-week Validation", "Full Time Series"))), size = 2, stat = "identity") +
#   geom_hline(yintercept = 0) + 
#   theme_bw() +
#   labs(x = "Date", y = "Percent Error [%]", color = "Model") +
#   scale_fill_manual(name = "Model", 
#                     values = c("2-week Validation" = brewer.pal(8,"Dark2")[6], "Full Time Series" = brewer.pal(8,"Dark2")[3]),
#                     labels = c("2-week Validation", "Full Time Series"),
#                     guide = "legend") +
#   theme(legend.position = "none") +
#   theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
#   theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold"))


ggplot(data = cross.val.long) +
  geom_histogram(aes(abs(my.perc.error*100), fill = factor(month(Date), levels = c(1:12))), color = "black", binwidth = 2) +
  theme_bw() +
  scale_fill_manual(values = scico(12, palette = "romaO", direction = -1))+
  labs(x = "Relative error (%)", y = "Frequency") +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  labs(tag = "B", fill = "Mean month of \ncross-validatoin \ndata subset") + theme(plot.tag = element_text(size = 18, face = "bold")) +
  annotate("point", x = 100*m3.rel.error, y = 3.5, color = "black", shape = 13, size = 3, stroke = 1) 


# FIGURE 4A
plot.4a <- ggplot(data = my.cross.val.values) +
  geom_point(aes(x = mean.O2bio, y = mean.O2biopred, color = month(Date)), size = 3) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.5, lwd = 1) +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white", color = NULL), strip.text = element_text(size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12)) +
  labs(tag = "A") + theme(plot.tag = element_text(size = 18, face = "bold")) +
  ylim(c(min(my.cross.val.values$mean.O2bio)*1.1, max(my.cross.val.values$mean.O2bio)*1.1)) + xlim(c(min(my.cross.val.values$mean.O2bio)*1.1, max(my.cross.val.values$mean.O2bio)*1.1)) +
  theme(panel.grid = element_blank()) +
  labs(color = "Mean month of \ncross-validatoin \ndata subset", x = "Mean [O2]bio (μM)", y = "Mean predicted [O2]bio (μM)") +
  scale_color_gradientn(colors = scico(12, palette = "romaO", direction = -1), labels = seq(0,12,3), breaks = seq(0,12,3), guide = guide_colorbar(reverse = T)) +
  annotate("point", x = mean(test$O2bio.estimated), y = mean(o2bio.predict$predictions), color = "black", shape = 13, size = 3, stroke = 1) +
  theme(legend.position = "none")
  
plot.4a

length(which(cross.val.long$my.perc.error <= m3.rel.error)) / nrow(cross.val.long)

test$O2bio.estimated - o2bio.predict$predictions













