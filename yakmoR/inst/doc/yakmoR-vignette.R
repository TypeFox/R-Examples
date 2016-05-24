## ---- echo=TRUE, results='asis', eval=TRUE-------------------------------
library(yakmoR)

data(iris)
irisM = as.matrix(iris[sample(nrow(iris)), -5]) # convert to matrix, also remove class-information
dat = irisM[1:100, ] # take first 100 data points for clustering
resObj = yakmoR::orthoKMeansTrain (x = dat, k = 3,  rounds = 4)
centers2 = resObj$centers[[2]] # centers of 2nd round

dat = as.matrix( irisM[101:nrow(irisM), -5]) # take rest of data for prediction
results = yakmoR::orthoKMeansPredict (x = dat, obj = resObj)

