## SIMCA missclassification

## Computing missclassification rates and confusion matrix for
##  SIMCA is available, but not very convenient and not documented.
##
##  If you call predict without 'newdata', resubstitution will be
##  performed and the apparent error rate and confusion table will be printed.
##
##  If you have new data and you know the group membership,
##  you can compute missclassification error and confusion
##  table as shown in the example below.
##

library(rrcovHD)
data(pottery)
dim(pottery)        # 27 observations in 2 classes, 6 variables
head(pottery)

## Build the SIMCA model. Use RSimca for a robust version
cs <- CSimca(origin~., data=pottery)
cs
summary(cs)

predict(cs)       ## without new data - perform resubstitution

## generate a sample from the pottery data set -
##  this will be the "new" data to be predicted
smpl <- sample(1:nrow(pottery), 5)
test <- pottery[smpl, -7]          # extract the test sample. Remove the last (grouping) variable
grp <- pottery[smpl, 7]            # remember the grouping variable also
print(test)


## predict new data
pr <- predict(cs, newdata=test)

pr@classification
ctab <- rrcov::mtxconfusion(grp, pr@classification)
ctab
acctab <- t(apply(ctab, 1, function(x) x/sum(x)))
dimnames(acctab) <- dimnames(ctab)
acctab

miss <- 1 - sum(diag(ctab))/sum(ctab)
miss
