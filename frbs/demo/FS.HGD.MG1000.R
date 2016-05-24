library(frbs)

## Input data
data(frbsData)
data.train <- frbsData$MackeyGlass1000.dt[1: 500, ]
data.fit <- data.train[, 1 : 4]
data.tst <- frbsData$MackeyGlass1000.dt[501 : 1000, 1 : 4]
real.val <- matrix(frbsData$MackeyGlass1000.dt[501 : 1000, 5], ncol = 1)
range.data<-matrix(c(0.43462, 1.3105, 0.43462, 1.3105, 0.43462, 1.3105, 0.43462, 1.3105, 0.43462, 1.3105), nrow=2)

## Set the method and its parameters
control <- list(num.labels = 5, max.iter = 100, step.size = 0.01, alpha.heuristic = 1, type.tnorm = "MIN", 
               type.snorm = "MAX", type.implication.func = "ZADEH", name = "MG1000")  
method.type <- "FS.HGD"

## Generate fuzzy model
object <- frbs.learn(data.train, range.data, method.type, control)

## Fitting step
res.fit <- predict(object, data.fit)

## Predicting step
res.test <- predict(object, data.tst)

## Error calculation
y.pred <- res.test
y.real <- real.val
bench <- cbind(y.pred, y.real)
colnames(bench) <- c("pred. val.", "real. val.")
print("Comparison FS.HGD Vs Real Value on Mackey Glass Data Set")
print(bench)

residuals <- (y.real - y.pred)
MSE <- mean(residuals^2)
RMSE <- sqrt(mean(residuals^2))
SMAPE <- mean(abs(residuals)/(abs(y.real) + abs(y.pred))/2)*100
err <- c(MSE, RMSE, SMAPE)
names(err) <- c("MSE", "RMSE", "SMAPE")
print("FS.HGD: Error Measurement: ")
print(err) 

## Comparing between simulation and real data
op <- par(mfrow = c(2, 1))
x1 <- seq(from = 1, to = nrow(res.fit))
result.fit <- cbind(data.train[, 5], res.fit)
plot(x1, result.fit[, 1], col="red", main = "Mackey Glass: Fitting phase (the training data(red) Vs Sim. result(blue))", type = "l", ylab = "MG")
lines(x1, result.fit[, 2], col="blue")

result.test <- cbind(real.val, res.test)
x2 <- seq(from = 1, to = nrow(result.test))
plot(x2, result.test[, 1], col="red", main = "Mackey Glass: Predicting phase (the Real Data(red) Vs Sim. result(blue))", type = "l", ylab = "MG")
lines(x2, result.test[, 2], col="blue", type = "l")

par(op)

