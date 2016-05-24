library(discSurv)
library(pec)
library(mgcv)

######################################
# tprUno, fprUno, aucUno, concorIndex

# Check if simple case of just training data without test sets works
# Example with unemployment data
data(cost)
summary(cost$time)
# Convert time to quarters
cost$time <- ceiling(cost$time/120)
summary(cost$time)
# Stratified sampling
cost1 <- subset(cost, status==1)
cost0 <- subset(cost, status==0)
set.seed(-567)
Ind1 <- sample(1:dim(cost1) [1], 50)
set.seed(-789)
Ind2 <- sample(1:dim(cost0) [1], 50)
# Combine
cost <- rbind(cost1 [Ind1, ], cost0 [Ind2, ])
head(cost)
summary(cost$status)
summary(cost$time)

# Estimate true positive rate of time interval 1: 
# Correspondes to two month duration (each interval is of length two weeks)
tryTPR <- tprUno (timepoint=1, dataSet=cost, trainIndices=list(1:dim(cost)[1]), survModelFormula=time ~ prevStroke + cholest, censModelFormula=status ~ 1, linkFunc="logit", idColumn=NULL, timeAsFactor=FALSE)
tryTPR
plot(tryTPR)

# Estimate false positive rate of time interval 1:
tryFPR <- fprUno (timepoint=1, dataSet=cost, trainIndices=list(1:dim(cost)[1]),  survModelFormula=time ~ prevStroke + cholest, censModelFormula=status ~ 1, linkFunc="logit", idColumn=NULL, timeAsFactor=FALSE)
tryFPR
plot(tryFPR)

# Estimate false positive rate of time interval 1:
tryAUC <- aucUno (tprObj=tryTPR, fprObj=tryFPR)
tryAUC
plot(tryAUC)

# Estimate global concordance index:
tryConcorIndex <- concorIndex (tryAUC)
tryConcorIndex
summary(tryConcorIndex)

################################################
# tprUnoShort, fprUnoShort:
# aucUno, concorIndex based on Short versions of tpr, fpr

set.seed(-367)
Ind1 <- sample(1:dim(cost1) [1], 100)
set.seed(-1089)
Ind2 <- sample(1:dim(cost0) [1], 100)
costTrain <- rbind(cost1 [Ind1 [1:50], ], cost0 [Ind2 [1:50], ])
costTest <- rbind(cost1 [Ind1 [51:100], ], cost0 [Ind2 [51:100], ])

# Train model
costTrainLong <- dataLong(dataSet=costTrain, timeColumn="time", censColumn="status")
costTestLong <- dataLong(dataSet=costTest, timeColumn="time", censColumn="status")
gamTrain <- gam(formula=y ~ timeInt + prevStroke + cholest, data=costTrainLong, family=binomial())
gamMarker <- predict(gamTrain, newdata=cbind(costTest, timeInt=factor(costTest$time)))

# Estimate true positive rate of time interval 1: 
# Correspondes to two month duration (each interval is of length two weeks)
tryTPRshort <- tprUnoShort (timepoint=1, marker=gamMarker, newTime=costTest$time, newEvent=costTest$status, trainTime=costTrain$time, trainEvent=costTrain$status)
tryTPRshort
plot(tryTPRshort)

# Estimate false positive rate of time interval 1:
tryFPRshort <- fprUnoShort (timepoint=1, marker=gamMarker, newTime=costTest$time, 
                            newEvent=costTest$status)
tryFPRshort
plot(tryFPRshort)

# Estimate false positive rate of time interval 1:
tryAUCshort <- aucUno (tprObj=tryTPRshort, fprObj=tryFPRshort)
tryAUCshort
plot(tryAUCshort)

# Estimate global concordance index:
tryConcorIndexShort <- concorIndex (tryAUCshort)
tryConcorIndexShort
summary(tryConcorIndexShort)

############################
# martingaleResid

# Calculate martingale residuals for the unemployment data subset
MartResid <- martingaleResid (dataSet=cost, survModelFormula=time ~ prevStroke + cholest, censColumn="status", linkFunc="logit", idColumn=NULL)
MartResid
CheckAlmostZero <- all.equal(sum(MartResid$Output$MartingaleResid), 0, tolerance = .Machine$double.eps ^ (1/3))
stopifnot(CheckAlmostZero)
stopifnot(length(MartResid$Output$MartingaleResid)==dim(cost) [1])

# Plot martingale residuals vs each covariate in the event interval
# Dotted line is a loess estimate
plot(MartResid)

##########################
# devResid

# Calculate deviance residuals for the unemployment data subset
devianceResidualsCheck <- devResid (dataSet=cost, survModelFormula=time ~ prevStroke + cholest, censColumn="status", linkFunc="logit", idColumn=NULL)
devianceResidualsCheck
stopifnot(length(devianceResidualsCheck$Output$DevResid)==dim(cost) [1])

###########################
# adjDevResid

adjDevianceResidualsCheck <- adjDevResid (dataSet=cost, survModelFormula=time ~ prevStroke + cholest, censColumn="status", linkFunc="logit", idColumn=NULL)
adjDevianceResidualsCheck
plot(adjDevianceResidualsCheck)
stopifnot(length(adjDevianceResidualsCheck$Output$AdjDevResid)==dim(cost) [1])

#####################################
# brierScore

tryBrierScore <- brierScore (dataSet=cost, trainIndices=list(1:dim(cost)[1]), survModelFormula=time ~ prevStroke + cholest, linkFunc="logit", censColumn="status", idColumn=NULL)
tryBrierScore
stopifnot(length(tryBrierScore)==dim(cost) [1])

#######################################
# predErrDisc

# # One time point
# tryPredErrDisc1 <- predErrDisc (timepoints=1, dataSet=cost, trainIndices=list(1:dim(cost)[1]), survModelFormula=time ~ prevStroke + cholest, censModelFormula=status ~ 1, linkFunc="logit", idColumn=NULL)
# tryPredErrDisc1
# summary(tryPredErrDisc1)
# # Multiple time points
# tryPredErrDisc2 <- predErrDisc (timepoints=3:5, dataSet=cost, trainIndices=list(1:dim(cost)[1]), survModelFormula=time ~ prevStroke + cholest, censModelFormula=status ~ 1, linkFunc="logit", idColumn=NULL)
# tryPredErrDisc2
# summary(tryPredErrDisc2)

#######################################
# predErrDiscShort

# Check if errors are zero in case the observed surv function is specified
InsertTimePoints <- 1:10
InsertTrainTime <- 1:10
InsertTrainEvent <- rep(1, 10)
InsertNewTime <- 1:10
InsertNewEvent <- rep(1, 10)
TheorySobs <- vector("list", length(InsertTimePoints))

# Order: Times -> Persons
for(i in 1:length(InsertTimePoints)) {
  TheorySobs [[i]] <- ifelse(InsertTimePoints [i] < InsertNewTime, 1, 0)
}

# Order: Persons -> Times
InsertEstSurv <- vector("list", length(InsertNewTime))
for (i in 1:length(InsertNewTime)) {
  InsertEstSurv [[i]] <- sapply(1:length(TheorySobs), function (x) TheorySobs [[x]] [i])
}

# Check
tryPredErrDisc1Short <- predErrDiscShort (timepoints=InsertTimePoints, estSurvList=InsertEstSurv, newTime=InsertNewTime, newEvent=InsertNewEvent, trainTime=InsertTrainTime, trainEvent=InsertTrainEvent)
stopifnot(all(tryPredErrDisc1Short$Output$predErr == 0))

# Case with known solutions
# Check if results are reproduceable
InsertTimePoints <- 1:100
InsertTrainTime <- 1:100
InsertTrainEvent <- rep(c(rep(0, 5), rep(1, 5)), 10)
InsertNewTime <- 30:49
InsertNewEvent <- rep(c(rep(1, 5), rep(0, 5)), 2)

# Order: Times -> Persons
TheorySobs <- vector("list", length(InsertTimePoints))
for(i in 1:length(InsertTimePoints)) {
  TheorySobs [[i]] <- ifelse(InsertTimePoints [i] < InsertNewTime, 1, 0)
}

# Order: Persons -> Times
TheorySobsOrd <- vector("list", length(InsertNewTime))
InsertEstSurv <- vector("list", length(InsertNewTime))
squDiff <- matrix(0, nrow=length(InsertNewTime), ncol=length(InsertTimePoints))
for (i in 1:length(InsertNewTime)) {
  TheorySobsOrd [[i]] <- sapply(1:length(TheorySobs), function (x) TheorySobs [[x]] [i])
  # Construct estimations
  InsertEstSurv [[i]] <- rep(0.5, length(TheorySobs))
  # Calculate squared differences
  squDiff [i, ] <- (TheorySobsOrd [[i]] - InsertEstSurv [[i]])^2
}

# Estimate survival function of censoring process
LongTrainData <- dataLong(dataSet=data.frame(InsertTrainTime=InsertTrainTime, InsertTrainEvent=InsertTrainEvent), timeColumn="InsertTrainTime", censColumn="InsertTrainEvent")
LongTrainData <- dataCensoring(dataSet=LongTrainData, respColumn="y", idColumn="obj")
LongTrainData <- na.omit(LongTrainData)
glmFit <- glm(yCens ~ timeInt, data=LongTrainData, family=binomial(), control=glm.control(maxit = 2500))
factorPrep <- factor(1:max(as.numeric(as.character(LongTrainData$timeInt))))
Gpre <- cumprod(1 - predict(glmFit, newdata=data.frame(timeInt=factorPrep), type="response"))
Gtheory <- c(1, Gpre)
# Iterate over all time points
results <- vector("numeric", length(InsertTimePoints))
Weights <- vector("list", length(InsertTimePoints))
for(i in 1:length(InsertTimePoints)) {
  Weights [[i]] <- as.numeric(InsertNewEvent * (1 - TheorySobs [[InsertTimePoints [i] ]]) / Gtheory [InsertNewTime]) + (TheorySobs [[InsertTimePoints [i] ]] / Gpre [InsertTimePoints [i] ])
  results [i] <- sum(squDiff [, i] * Weights [[i]]) / length(Weights [[i]])
}
names(results) <- paste("T=", InsertTimePoints, sep="")
RemoveIndex <- is.nan(results) | is.na(results) | is.infinite (results)
results <- results [!RemoveIndex]

# Check
tryPredErrDisc1Short <- predErrDiscShort (timepoints=InsertTimePoints, estSurvList=InsertEstSurv, newTime=InsertNewTime, newEvent=InsertNewEvent, trainTime=InsertTrainTime, trainEvent=InsertTrainEvent)
stopifnot(all.equal(sum(abs(results - tryPredErrDisc1Short$Output$predErr)), 0))

######################################
# intPredErrDisc

# # Long version
# checkIntPredErrDisc <- intPredErrDisc (tryPredErrDisc1)
# checkIntPredErrDisc
# stopifnot(checkIntPredErrDisc>=0)

# Short version
checkIntPredErrDisc2 <- intPredErrDisc (predErrObj=tryPredErrDisc1Short)
checkIntPredErrDisc2
stopifnot(checkIntPredErrDisc2>=0)

# Short version up to time point 3
checkIntPredErrDisc2 <- intPredErrDisc (predErrObj=tryPredErrDisc1Short, tmax=3)
checkIntPredErrDisc2
stopifnot(checkIntPredErrDisc2>=0)
