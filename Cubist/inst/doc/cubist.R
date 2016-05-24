### R code from vignette source 'cubist.Rnw'

###################################################
### code chunk number 1: startup
###################################################
library(mlbench)
data(BostonHousing)
library(Cubist)


###################################################
### code chunk number 2: bh1
###################################################
library(Cubist)
library(mlbench)
data(BostonHousing)
BostonHousing$chas <- as.numeric(BostonHousing$chas) - 1

set.seed(1)

inTrain <- sample(1:nrow(BostonHousing), floor(.8*nrow(BostonHousing)))

trainingPredictors <- BostonHousing[ inTrain, -14]
testPredictors     <- BostonHousing[-inTrain, -14]

trainingOutcome <- BostonHousing$medv[ inTrain]
testOutcome     <- BostonHousing$medv[-inTrain]

modelTree <- cubist(x = trainingPredictors, y = trainingOutcome)
modelTree


###################################################
### code chunk number 3: bh2
###################################################
summary(modelTree)


###################################################
### code chunk number 4: bh3
###################################################
mtPred <- predict(modelTree, testPredictors)
## Test set RMSE
sqrt(mean((mtPred - testOutcome)^2))
## Test set R^2
cor(mtPred, testOutcome)^2


###################################################
### code chunk number 5: bh4
###################################################
set.seed(1)
committeeModel <- cubist(x = trainingPredictors, y = trainingOutcome,
                         committees = 5)
summary(committeeModel)


###################################################
### code chunk number 6: bh5
###################################################
cmPred <- predict(committeeModel, testPredictors)
## RMSE
sqrt(mean((cmPred - testOutcome)^2))
## R^2
cor(cmPred, testOutcome)^2


###################################################
### code chunk number 7: bh6
###################################################
instancePred <- predict(committeeModel, testPredictors, neighbors = 5)
## RMSE
sqrt(mean((instancePred - testOutcome)^2))
## R^2
cor(instancePred, testOutcome)^2


###################################################
### code chunk number 8: tune
###################################################
library(caret)

set.seed(1)
cTune <- train(x = trainingPredictors, y = trainingOutcome,
               "cubist",
               tuneGrid = expand.grid(.committees = c(1, 10, 50, 100), 
                                      .neighbors = c(0, 1, 5, 9)),
               trControl = trainControl(method = "cv"))
cTune


###################################################
### code chunk number 9: tune
###################################################
trellis.par.set(caretTheme())
print(plot(cTune, aut.key = list(columns = 4)))


###################################################
### code chunk number 10: lstat
###################################################
lstat <- trainingPredictors[, "lstat", drop = FALSE]
justRules <- cubist(lstat, trainingOutcome)
andCommittees <- cubist(lstat, trainingOutcome, committees = 100)


###################################################
### code chunk number 11: lstatPlot
###################################################
lstatTest <- testPredictors[, "lstat", drop = FALSE]
newOrder <- order(lstatTest$lstat)
lstatTest <- lstatTest[newOrder,,drop = FALSE]
testOutcome <- testOutcome[newOrder]
plot(lstatTest$lstat, testOutcome,
     pch = 16, col = rgb(.2, .2, .2, .5),
     xlab = "lstat", ylab = "Median Home Value")
points(lstatTest$lstat, predict(justRules, lstatTest),
       type = "l", lwd = 2, col = "black")
points(lstatTest$lstat, predict(justRules, lstatTest, neighbors = 5),
       type = "l", lwd = 2, col = "blue")
points(lstatTest$lstat, predict(andCommittees, lstatTest),
       type = "l", lwd = 2, col = "darkred")
legend(20, 50,
       c("Rules", "100 Committees", "Rules + 5 Neighbors"),
       col = c("black", "darkred", "blue"),
       lwd = rep(2, 3))


###################################################
### code chunk number 12: vimp
###################################################
summary(modelTree)
modelTree$usage
library(caret)
varImp(modelTree)


###################################################
### code chunk number 13: <session
###################################################
toLatex(sessionInfo())


