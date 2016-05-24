##  SIMCA is basically a set of several PCA models (as many as the classes).
##
##  1. getLoadings is a function of the PCA object
##  2. If you extract the PCA objects of the CSimca object,
##      you can do all its plots. To do this use obj@pcaobj[[i]]
##      where obj is a CSimca or RSimca object and i is the number of the class.
##

######################################################
library(rrcovHD)
data(pottery)
dim(pottery)        # 27 observations in 2 classes, 6 variables
head(pottery)

## Build the SIMCA model. Use RSimca for a robust version
cs <- CSimca(origin~., data=pottery)
cs
summary(cs)

str(cs)     # you will se the element pcaobj:List of 2

## With these PCA objects you can do whatever you can do with any other PCA object
cs@pcaobj[[1]]           # the PCA corresponding to the first class
cs@pcaobj[[2]]           # the PCA corresponding to the second class

## Extract the loadings of the first and second PCA object
## Unfortunately there is no loadings plot yet, but you could program it yourself, having the loadings
getLoadings(cs@pcaobj[[1]])
getLoadings(cs@pcaobj[[2]])

## PCA plots
screeplot(cs@pcaobj[[2]])
biplot(cs@pcaobj[[2]])
plot(cs@pcaobj[[2]])
