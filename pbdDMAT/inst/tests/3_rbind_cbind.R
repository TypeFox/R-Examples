library(pbdTEST)
settings(mpi=TRUE)

.BLDIM <- 2
comm.set.seed(seed=1234, diff=FALSE)


### --------------------------------------
module("r/c-bind")

n <- 1e2
p <- 25

x <- matrix(rnorm(n*p), n, p)
y <- matrix(rnorm(n*p, mean=100, sd=10), n, p)

dx <- as.ddmatrix(x)
dy <- as.ddmatrix(y)

test("rbind()", {
  a <- rbind(x, y)
  b <- as.matrix(rbind(dx, dy))
})

test("cbind()", {
  a <- cbind(x, y)
  b <- as.matrix(cbind(dx, dy))
})

collect()



finalize()

