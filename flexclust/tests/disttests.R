library(flexclust)

## create a random matrix
x <- matrix(rnorm(1000), ncol=4)
rownames(x) <- 1:nrow(x)

## test canberra at margin
x[1:10,1] <- 0
x[1:20,2] <- 0
x[1:30,3] <- 0

stopifnot(all.equal(dist2(x,x,"eucl"), as.matrix(dist(x, "eucl"))))
stopifnot(all.equal(dist2(x,x,"man"), as.matrix(dist(x, "man"))))
stopifnot(all.equal(dist2(x,x,"max"), as.matrix(dist(x, "max"))))
stopifnot(all.equal(dist2(x,x,"bin"), as.matrix(dist(x, "bin"))))
stopifnot(all.equal(dist2(x,x,"can"), as.matrix(dist(x, "can"))))

stopifnot(all.equal(dist2(x,x,"eucl"), dist2(x,x,"mink",p=2)))
stopifnot(all.equal(dist2(x,x,"man"), dist2(x,x,"mink",p=1)))
