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
method.type <- "GFS.GCCL" 
control <- list(popu.size = 30, num.class = 3, num.labels = 3, persen_cross = 0.9, 
                     max.gen = 200, persen_mutant = 0.3,
                     name="sim-0") 

## Generate fuzzy model
object <- frbs.learn(tra.iris, range.data, method.type, control)

## Write frbsPMML file
## In this step, we provide two ways as follows.
## a. by calling frbsPMML() function directly. 
## b. by calling write.frbsPMML() function. 
####################

## a. by calling frbsPMML(), the result will be displayed in R console
objPMML <- frbsPMML(object)
 
## b. by calling write.frbsPMML(), the result will be saved as a file
##     in the working directory.
write.frbsPMML(objPMML, file = "GCCL.Iris")

## Read frbsPMML file
##############################
 
object.pmml <- read.frbsPMML("GCCL.Iris.frbsPMML")
 
## Perform predicting step
###############################

## Predicting step
res.test <- predict(object.pmml, tst.iris)

## error calculation
err = 100*sum(real.iris!=res.test)/nrow(real.iris)

print("The result: ")
print(res.test)
print("GFS-GCCL: percentage Error on Iris: ")
print(err) 



