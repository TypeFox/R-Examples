library(pbdTEST)
settings(mpi=TRUE)

.BLDIM <- 2
comm.set.seed(seed=1234, diff=FALSE)


### --------------------------------------
module("ddmatrix constructor")

n <- 1e2
p <- 25


y <- matrix(rnorm(n*p, mean=100, sd=10), n, p)
### Need to do y first, since otherwise the rng's separate...
if (comm.rank()==0){
  x <- matrix(rnorm(n*p), n, p)
} else {
  x <- NULL
}

dx <- as.ddmatrix(x)
dy <- as.ddmatrix(y)

test("as.matrix --- all ranks", {
  a <- x
  b <- as.matrix(dx)
})

test("as.matrix --- proc.dest=0", {
  a <- y
  b <- as.matrix(dy, proc.dest=0)
})

collect()



finalize()

