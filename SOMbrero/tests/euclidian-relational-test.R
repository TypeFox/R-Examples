# Check the equivalence between numeric SOM and  relational SOM for a 
# four-dimensional numeric dataset and the square of Euclidean distance used
# as dissimilarity

library(SOMbrero)

set.seed(123)
nsom <- trainSOM(iris[1:30,1:4], init.proto="obs", scaling="none", maxit=50)

set.seed(123)
iris.dist <- dist(iris[1:30,1:4], method="euclidian", diag=TRUE, upper=TRUE)^2
rsom <- trainSOM(x.data=iris.dist, type="relational", scaling="none", maxit=50)

stopifnot(identical(nsom$clustering, rsom$clustering))