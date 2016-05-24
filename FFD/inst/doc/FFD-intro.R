### R code from vignette source 'FFD-intro.Rnw'

###################################################
### code chunk number 1: FFD-intro.Rnw:31-37
###################################################
set.seed(1504)
#options(width=70, prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)
options(width=70, prompt = "> ", continue = "   ", useFancyQuotes = FALSE)
ps.options(family="Times")
library("FFD")
data("sheepData")


###################################################
### code chunk number 2: FFD-intro.Rnw:161-163
###################################################
nVec <- 1:300
outVec <- computeAlpha(nAnimalVec = nVec, method = "individual", herdSensitivity = 0.7, intraHerdPrevalence = 0.2, diagSensitivity = 0.9)


###################################################
### code chunk number 3: FFD-intro.Rnw:169-171
###################################################
plot(nVec, 1-outVec, type = "l", xlab = "Herd size", ylab = "Herd sensitivity")
abline(h = 0.7, lty = 2, lwd = 2, col = "red")


###################################################
### code chunk number 4: FFD-intro.Rnw:563-569
###################################################
data(sheepData)
mySurvey <- surveyData(nAnimalVec = sheepData$nSheep,
    populationData = sheepData, designPrevalence = 0.002,
    alpha = 0.05, intraHerdPrevalence = 0.2,
    diagSensitivity = 0.9, costHerd = 30, costAnimal = 7)
summary(mySurvey)


###################################################
### code chunk number 5: FFD-intro.Rnw:584-587
###################################################
myIndSamplingSummary <- indSamplingSummary(survey.Data = mySurvey,
    stepSize = 0.05)
summary(myIndSamplingSummary)


###################################################
### code chunk number 6: FFD-intro.Rnw:600-601
###################################################
plot(myIndSamplingSummary)


###################################################
### code chunk number 7: FFD-intro.Rnw:621-624
###################################################
myIndSampling <- indSampling(survey.Data = mySurvey,
    herdSensitivity = 0.7)
summary(myIndSampling)


###################################################
### code chunk number 8: FFD-intro.Rnw:637-640
###################################################
myLtdSampleSummary <- ltdSamplingSummary(survey.Data = mySurvey,
    sampleSizeLtdMax = 30)
summary(myLtdSampleSummary)


###################################################
### code chunk number 9: FFD-intro.Rnw:653-654
###################################################
plot(myLtdSampleSummary)


###################################################
### code chunk number 10: FFD-intro.Rnw:674-676
###################################################
myLtdSampling <- ltdSampling(survey.Data = mySurvey, sampleSizeLtd = 7)
summary(myLtdSampling)


###################################################
### code chunk number 11: FFD-intro.Rnw:685-704
###################################################
## Fixed sampling:
##################
sampleFixed <- sample(x = myIndSampling, size = "fixed")
## Sample Size:
length(sampleFixed$indexSample)
## Significance:
sampleFixed$aPostAlpha
## Sample:
head(sampleFixed$sample)

## Dynamic sampling:
####################
sampleDynamic <- sample(x = myIndSampling, size = "dynamic")
## Sample Size:
length(sampleDynamic$indexSample)
## Significance:
sampleDynamic$aPostAlpha
## Sample:
head(sampleFixed$sample)


###################################################
### code chunk number 12: FFD-intro.Rnw:713-717
###################################################
p.value <- computePValue(nPopulation = 15287, nSample = 1630,
    nDiseased = round(15287*0.002), sensitivity = 0.8633,
    specificity = 1)
p.value


###################################################
### code chunk number 13: FFD-intro.Rnw:722-726
###################################################
nSample <- computeOptimalSampleSize(nPopulation = 15287,
    prevalence = 0.002, alpha = 0.05, sensitivity = 0.8633,
    specificity = 1, lookupTable = FALSE)
nSample


###################################################
### code chunk number 14: FFD-intro.Rnw:737-741
###################################################
lookupTable <- computeOptimalSampleSize(nPopulation = max(sheepData$nSheep),
    prevalence = 0.2, alpha = 0.3, sensitivity = 0.9, specificity = 1,
    lookupTable = TRUE)
lookupTable


###################################################
### code chunk number 15: FFD-intro.Rnw:751-756
###################################################
alphaList <- computeAlphaLimitedSampling(stockSizeVector = sheepData$nSheep,
    sampleSizeLtd = 7, intraHerdPrevalence = 0.2, diagSensitivity = 0.9,
    diagSpecificity = 1)
str(alphaList$alphaDataFrame)
alphaList$meanAlpha


###################################################
### code chunk number 16: FFD-intro.Rnw:769-784
###################################################
sampleVec <- sample(sheepData$nSheep, 2550, replace = FALSE)
alphaVec <- computeAlpha(nAnimalVec = sampleVec, method = "limited",
    sampleSizeLtd = 9, intraHerdPrevalence = 0.2, diagSensitivity = 0.9)

system.time({
errorExact <- computeAposterioriError(alphaErrorVector = alphaVec,
    nPopulation = 5000, nDiseased = 5, method = "exact")
})
errorExact

system.time({
errorApprox <- computeAposterioriError(alphaErrorVector = alphaVec,
    nPopulation = 5000, nDiseased = 5, method = "approx")
})
errorApprox


