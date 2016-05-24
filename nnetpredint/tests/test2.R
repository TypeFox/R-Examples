# Example 2: Using the nnet object trained by nnet package
library(nnet)
xTrain <- rbind(cbind(runif(150,min = 0, max = 0.5),runif(150,min = 0, max = 0.5)) ,
		cbind(runif(150,min = 0.5, max = 1),runif(150,min = 0.5, max = 1))
		)
nObs <- dim(xTrain)[1]
yTrain <- 0.5 + 0.4 * sin(2* pi * xTrain %*% c(0.4,0.6)) +rnorm(nObs,mean = 0, sd = 0.05)
plot(xTrain %*% c(0.4,0.6),yTrain)

# Training nnet models
net <- nnet(yTrain ~ xTrain,size = 3, rang = 0.1,decay = 5e-4, maxit = 500)
yFit <- c(net$fitted.values)
nodeNum <- c(2,3,1)
wts <- net$wts

# New data for prediction intervals
library(nnetpredint)
newData <- cbind(seq(0,1,0.05),seq(0,1,0.05))
yTest <- 0.5 + 0.4 * sin(2* pi * newData %*% c(0.4,0.6))+rnorm(dim(newData)[1],mean = 0, sd = 0.05)

# S3 generic method: Object of nnet
yPredInt <- nnetPredInt(net, xTrain, yTrain, newData)
print(yPredInt[1:20,])

# S3 default method for user defined input
yPredInt2 <- nnetPredInt(object = NULL, xTrain, yTrain, yFit, node = nodeNum, wts = wts, newData, 
	alpha = 0.05, funName = 'sigmoid')

plot(newData %*% c(0.4,0.6),yTest,type = 'b')
lines(newData %*% c(0.4,0.6),yPredInt$yPredValue,type = 'b',col='blue')
lines(newData %*% c(0.4,0.6),yPredInt$lowerBound,type = 'b',col='red')   # lower bound
lines(newData %*% c(0.4,0.6),yPredInt$upperBound,type = 'b',col='red')   # upper bound