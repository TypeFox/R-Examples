##  SIMCA example
##  Using CSimca (non robust) for classifying new observations.
##

library(rrcovHD)
library(plsgenomics)
data(Colon)

names(Colon)

## Colon$Y is the grouping variable (1=Normal, 2=Tumor) and
## Colon$X contains the data, 62x2000, 62 samples in 2000 variables
##

## Transform the data (see Pires and Branco, 2010)
##
data1 <- log(Colon$X)
med.row <- apply(data1, 1, median)
mad.row <- apply(data1, 1, mad)
data2 <- sweep(data.matrix(data1), 1, med.row)    # subtract the row medians
data2 <- sweep(data2, 1, mad.row, FUN="/")        # devide by the row MADS

colon <- cbind.data.frame(data2, grp=Colon$Y)

## Build the models (one PCA model for each group)
## Choose 8 components in each group (k=8)

sim <- RSimca(grp~., data=colon, k=8)
sim

## Use the 'sim' object to predict new observations
##  for the sake of the example take the first 10
##  observations from the training data set

test <- as.matrix(colon[1:10,-2001])
predict(sim, newdata=test)
