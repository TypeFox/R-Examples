stopifnot(require(hdi))

## this is the example code of the help file of plotClusterGroupBound


## Create a regression problem with correlated design (n = 10, p = 3):
## a block of size 2 and a block of size 1, within-block correlation is 0.99

set.seed(29)
p   <- 3
n   <- 10

Sigma <- diag(p)
Sigma[1,2] <- Sigma[2,1] <- 0.99

x <- matrix(rnorm(n * p), nrow = n) %*% chol(Sigma)

## Create response with active variable 1
beta    <- rep(0, p)
beta[1] <- 5

y  <- as.numeric(x %*% beta + rnorm(n))

## Compute the lower bound for all groups in a hierarchical clustering tree
cgb5 <- clusterGroupBound(x, y, nsplit = 4) ## use larger value for nsplit!

## Plot the tree with y-axis proportional to the (log) of the number of
## group members and node sizes proportional to the lower l1-norm bound.
plot(cgb5)

## Show the lower bound on the y-axis and node sizes proportional to
## number of group members
plot(cgb5, yaxis = "")
