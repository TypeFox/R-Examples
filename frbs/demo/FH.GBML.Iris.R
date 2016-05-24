library(frbs)

## Input data
data(iris)
set.seed(2)
irisShuffled <- iris[sample(nrow(iris)),]
irisShuffled[,5] <- unclass(irisShuffled[,5])
tra.iris <- irisShuffled[1:105,]
tst.iris <- irisShuffled[106:nrow(irisShuffled),1:4]
real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)
range.data <- matrix(c(4.3, 7.9, 2.0, 4.4, 1.0, 6.9, 0.1, 2.5), nrow=2)

## Set the method and its parameters
method.type <- "FH.GBML" 
control <- list(popu.size = 10, max.num.rule = 50, num.class = 3, 
			persen_cross = 0.9, max.gen = 20, persen_mutant =0.3, p.dcare = 0.5, p.gccl = 0.5, 
			name="sim-0") 
 
## Generate fuzzy model
object <- frbs.learn(tra.iris, range.data, method.type, control)

## Predicting step
res.test <- predict(object, tst.iris)

## error calculation
err = 100*sum(real.iris!=res.test)/nrow(real.iris)

print("The result: ")
print(res.test)
print("FH-GBML: percentage Error on Iris: ")
print(err) 



