library(frbs)

## input data
data(frbsData)
data.train <- frbsData$GasFurnance.dt[1 : 204, ]
data.fit <- data.train[, 1 : 2]
data.tst <- frbsData$GasFurnance.dt[205 : 292, 1 : 2]
real.val <- matrix(frbsData$GasFurnance.dt[205 : 292, 3], ncol = 1)

range.data<-matrix(c(-2.716, 2.834, 45.6, 60.5, 45.6, 60.5), nrow=2)

## set the method and its parameters
method.type <- "ANFIS"
control <- list(num.labels = 3, max.iter = 100, step.size = 0.01, type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH", name = "GasFur")

## generate fuzzy model
object <- frbs.learn(data.train, range.data, method.type, control)

## Write frbsPMML file
## In this step, we provide two ways as follows.
## a. by calling frbsPMML() function directly. 
## b. by calling write.frbsPMML() function. 
####################

## a. by calling frbsPMML(), the result will be displayed in R console
objPMML <- frbsPMML(object)
 
## b. by calling write.frbsPMML(), the result will be saved as a file
##     in the working directory.
write.frbsPMML(objPMML, file = "ANFIS.GasFur")

## Read frbsPMML file
##############################
 
object.pmml <- read.frbsPMML("ANFIS.GasFur.frbsPMML")

## This process is a part of fitting the model using data training. 
res.fit <- predict(object.pmml, data.fit)

## predicting step
res.test <- predict(object.pmml, data.tst)

## error calculation
y.pred <- res.test
y.real <- real.val
bench <- cbind(y.pred, y.real)
colnames(bench) <- c("pred. val.", "real. val.")
print("Comparison ANFIS Vs Real Value on Gas Furnace Data Set")
print(bench)

residuals <- (y.real - y.pred)
MSE <- mean(residuals^2)
RMSE <- sqrt(mean(residuals^2))
SMAPE <- mean(abs(residuals)/(abs(y.real) + abs(y.pred))/2)*100
err <- c(MSE, RMSE, SMAPE)
names(err) <- c("MSE", "RMSE", "SMAPE")
print("Error Measurement: ")
print(err) 

## Comparing between simulation and real data
op <- par(mfrow = c(2, 1))
result.fit <- cbind(data.train[, 3], res.fit)
x1 <- seq(from = 1, to = nrow(result.fit))
plot(x1, result.fit[, 1], col="red", main = "Gas Furnance: Fitting phase (the training data(red) Vs Sim. result(blue))", type = "l", ylab = "CO2")
lines(x1, result.fit[, 2], col="blue", type = "l")

result.test <- cbind(real.val, res.test)
x2 <- seq(from = 1, to = nrow(result.test))
plot(x2, result.test[, 1], col="red", main = "Gas Furnance: Predicting phase (the Real Data(red) Vs Sim. result(blue))", type = "l", ylab = "CO2")
lines(x2, result.test[, 2], col="blue", type = "l")
par(op)

