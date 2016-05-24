##
##  SIMCA prediction
##
##  In order to predict with model build by CSimca
##  (or RSimca for the robust version) you can use
##  the standard function 'predict()'.
##  See the example below

library(rrcoHD)
data(pottery)
dim(pottery)        # 27 observations in 2 classes, 6 variables
head(pottery)

## Build the SIMCA model. Use RSimca for a robust version
cs <- CSimca(origin~., data=pottery)
cs
summary(cs)


## generate a sample from the pottery data set -
##  this will be the "new" data to be predicted
smpl <- sample(1:nrow(pottery), 5)
test <- pottery[smpl, -7]          # extract the test sample. Remove the last (grouping) variable
print(test)


## predict new data
pr <- predict(cs, newdata=test)

pr@classification 
