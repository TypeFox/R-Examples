library(frbs)

## Input data
data(iris)
set.seed(2)
irisShuffled <- iris[sample(nrow(iris)),]
irisShuffled[,5] <- unclass(irisShuffled[,5])
tra.iris <- irisShuffled[1:105,]
tst.iris <- irisShuffled[106:nrow(irisShuffled),1:4]
real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)
range.data.input <- matrix(c(4.3, 7.9, 2.0, 4.4, 1.0, 6.9, 0.1, 2.5), nrow=2)
 
## Set the method and its parameters
method.type <- "FRBCS.CHI"
control <- list(num.labels = 6, type.mf = "GAUSSIAN", type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH")

## Generate fuzzy model
object <- frbs.learn(tra.iris, range.data.input, method.type, control)

## Predicting step
res.test <- predict(object, tst.iris)

## error calculation
err = 100*sum(real.iris!=res.test)/nrow(real.iris)

print("The result: ")
print(res.test)

print("FRBCS.CHI: percentage Error on Iris")
print(err) 

