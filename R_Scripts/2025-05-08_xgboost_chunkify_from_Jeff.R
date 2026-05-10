set.seed(123)

library(caret)

combo <- readRDS('2026-05-06_combo.rds')

combo <- na.omit(combo)

combo$Date <- NULL

#### variable selection ####

## variable selection significantly degrades performance

# bad.predictors <- c(nearZeroVar(combo))
# bad.predictors <- c(bad.predictors, grep('lingulo', colnames(combo)))
# bad.predictors <- c(bad.predictors, grep('closterium', colnames(combo)))
# colnames(combo)[bad.predictors]
# combo[,bad.predictors] <- NULL

plot(combo$O2bio.estimated~combo$O2bio.predicted)

combo$O2bio.estimated <- NULL
combo$O2bio.predicted.altered <- NULL
response <- "O2bio.predicted"

#### transformation ####

## transformation has minimal (log10) or negative (RA) impact on performance

library(vegan)

#cols_to_transform <- setdiff(colnames(combo), response)

## log10

#combo[, cols_to_transform] <- log10(combo[, cols_to_transform])

## relative abundance

#temp <- decostand(combo[, cols_to_transform], 'total')
#combo <- cbind(temp, combo[,response])
#colnames(combo)[length(colnames(combo))] <- response 

## eliminate any inf

# combo[sapply(combo, is.numeric)] <- lapply(
#   combo[sapply(combo, is.numeric)], 
#   function(x) replace(x, is.infinite(x), 0)
# )

#### take a look at variables ####

library(gplots)

heatmap.2(data.matrix(combo[,colnames(combo) != response]),
          Rowv = F,
          dendrogram = 'none',
          margins = c(15,5))

#### model it ####

library(xgboost)

n_iter <- dim(combo)[1] - 20
all.cor <- vector(length = n_iter)
all.predictions <- vector("list", n_iter)

for(i in 1:(dim(combo)[1] - 20)){
  
  test <- combo[i:(i+19),]
  train <- combo[-c(i:(i+19)),]
  
  temp.xgb <- xgboost(
    data = train[,colnames(train)[grep(response, colnames(train), invert = T)]],
    y = train[[response]],
    nrounds = 100,
    verbose = 0
  )
  
  temp.prediction <- predict(temp.xgb, test)
  
  all.predictions[[i]] <- data.frame(
    predicted = temp.prediction,
    observed  = test[[response]],
    iteration = i
  )
  
  plot(temp.prediction, test$O2bio.predicted)
  print(cor(temp.prediction, test$O2bio.predicted))
  all.cor[i] <- cor(temp.prediction, test$O2bio.predicted)
}

mean(all.cor)

pred.df <- do.call(rbind, all.predictions)
plot(pred.df$predicted ~ pred.df$observed)
plot(all.cor,
     type = 'l')




