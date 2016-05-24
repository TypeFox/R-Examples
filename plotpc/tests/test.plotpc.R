# test.plotpc.R
# This is just a dummy test. Full tests can be found in tests/slowtests.
library(grid)
library(plotpc)
data(iris)
x <- iris[,c(3,4)] # select Petal.Length and Petal.Width
plotpc(x)
