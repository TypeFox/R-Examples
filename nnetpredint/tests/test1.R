# Example 1: Using the nn object trained by neuralnet package
set.seed(500)
library(MASS)
data <- Boston
maxs <- apply(data, 2, max)
mins <- apply(data, 2, min)
scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins)) # normalization
index <- sample(1:nrow(data),round(0.75*nrow(data)))
train_ <- scaled[index,]
test_ <- scaled[-index,]

library(neuralnet) # Training
n <- names(train_)
f <- as.formula(paste("medv ~", paste(n[!n %in% "medv"], collapse = " + ")))
nn <- neuralnet(f,data = train_,hidden = c(5,3),linear.output = FALSE)
plot(nn)

library(nnetpredint) # Getting Prediction confidence interval
x <- train_[,-14]
y <- train_[,14]
newData <- test_[,-14]

# S3 generic method: Object of nn
yPredInt <- nnetPredInt(nn, x, y, newData)
print(yPredInt[1:20,])

# S3 default method for user defined weight input
yFit <- c(nn$net.result[[1]])
nodeNum <- c(13,5,3,1)
m <- 3
wtsList <- nn$weights[[1]]
wts <- transWeightListToVect(wtsList,m)
yPredInt2 <- nnetPredInt(object = NULL, x, y, yFit, nodeNum, wts, newData, alpha = 0.05)
print(yPredInt2[1:20,])

# Compare to the predict values from the neuralnet Compute method
predValue <- compute(nn,newData)
print(matrix(predValue$net.result[1:20]))
