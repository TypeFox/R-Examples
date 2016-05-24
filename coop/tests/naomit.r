library(coop)
naomit <- coop:::naomit

set.seed(1234)

m <- 1000
n <- 20
len <- m*n

x <- matrix(rnorm(m*n, sd=10000), m, n)

y <- x
prop <- .01
y[sample(len, size=len*prop)] <- NA
# y[m,2] = NA

if (any(dim(na.omit(y)) == 0)) stop("zeros")
stopifnot(all.equal(na.omit(y), naomit(y), check.attributes=FALSE))

storage.mode(y) <- "integer"
stopifnot(all.equal(na.omit(y), naomit(y), check.attributes=FALSE))




m <- 10
n <- 2
len <- m*n

x <- matrix(rnorm(m*n, sd=10000), m, n)

y <- x
prop <- .01
y[sample(len, size=len*prop)] <- NA

if (any(dim(na.omit(y)) == 0)) stop("zeros")
stopifnot(all.equal(na.omit(y), naomit(y), check.attributes=FALSE))

storage.mode(y) <- "integer"
stopifnot(all.equal(na.omit(y), naomit(y), check.attributes=FALSE))



### TODO
# if (require(slam))
# {
#   library(slam)
#   csc <- as.simple_triplet_matrix(y)
#   z <- as.matrix(coop:::naomit_coo(as.double(csc$v), csc$i, csc$j))
#   stopifnot(all.equal(na.omit(y), z, check.attributes=FALSE))
# }
