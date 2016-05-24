## recreating exact fit problem for categorical data
require(robustbase)

## some simple balanced dataset with one grouping variable
ngrp <- 10
nrep <- 10
set.seed(2)
data <- data.frame(y = rnorm(ngrp*nrep), grp=rep(letters[1:ngrp], each=nrep))
## this works fine
m1 <- lmrob(y ~ grp, data)

## now contaminate the dataset
data2 <- data
data2$y[1:48] <- 1e10
try(m2 <- lmrob(y ~ grp, data2, trace.lev = 3))
## All observations of group "e" get rob. weight of 0:
weights <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) ## from trace output
weights %*% m1$x
