library(discSurv)

###########################################
# Test
# devResidShort (dataSet, hazards)
# adjDevResidShort (dataSet, hazards)

# Preprocessing
library(survival)
set.seed(0)
Indizes <- sample(1:dim(transplant)[1], 100)
transplant2 <- transplant[Indizes, ]
censColumns <- data.frame(model.matrix(~event, transplant2)[,-1])
transplant2 <- data.frame(transplant2, death=censColumns[,1])
transplant2$futime <- transplant2$futime+1
transplantLong <- dataLong(dataSet=transplant2, timeColumn="futime", censColumn="death")
transplantLong$timeInt <- as.numeric(as.character(transplantLong$timeInt))

# Estimate discrete hazards with gam
library(mgcv)
fitGAM <- gam(y ~ s(timeInt) + sex + abo + year, data=transplantLong, family=binomial())
hazPreds <- predict(fitGAM, type="response")

# Calculate deviance residuals
results <- devResidShort(dataSet=transplantLong, hazards=hazPreds)
results
results <- adjDevResidShort(dataSet=transplantLong, hazards=hazPreds)
results
