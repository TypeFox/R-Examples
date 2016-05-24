library(SOMbrero)

iris.dist <- dist(iris[1:30,1:4], method="minkowski", diag=TRUE, upper=TRUE, 
                  p=4)
rsom <- trainSOM(x.data=iris.dist, type="relational")

stopifnot(all.equal(as.vector(rowSums(rsom$prototypes)), 
                    rep(1, prod(rsom$parameters$the.grid$dim))))

stopifnot(!sum(rsom$prototypes<0))