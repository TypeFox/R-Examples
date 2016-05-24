library(pbdTEST)
settings(mpi=TRUE)

.BLDIM <- 2
comm.set.seed(seed=1234, diff=FALSE)

tol <- 1e-8

na_maker <- function(n)
{
  sapply(X=1:n, FUN=function(.) {
    cf <- sample(1:2, size=1, prob=c(.1, .9))
    if (cf==1) return(NA)
    else return(rnorm(1))
  })
}


### --------------------------------------
module("NA removal:  square")

for (i in 1:3){
  x <- matrix(na_maker(100), 10)
  dx <- as.ddmatrix(x)
  
  test(paste0("Random test #", i), {
    a <- na.exclude(x)
    b <- as.matrix(na.exclude(dx))
  }, check.attributes=FALSE)
}

collect()


### --------------------------------------
module("NA removal:  row")

for (i in 1:3){
  x <- matrix(na_maker(100), 1)
  dx <- as.ddmatrix(x)
  
  test(paste0("Random test #", i), {
    a <- na.exclude(x)
    b <- as.matrix(na.exclude(dx))
  }, check.attributes=FALSE)
}

for (i in 1:3){
  x <- matrix(na_maker(100), 1)
  dx <- as.ddmatrix(x, bldim=100)
  
  test(paste0("Random test (big bldim) #", i), {
    a <- na.exclude(x)
    b <- as.matrix(na.exclude(dx))
  }, check.attributes=FALSE)
}

collect()




### --------------------------------------
module("NA removal:  column")

for (i in 1:3){
  x <- matrix(na_maker(100), 100)
  dx <- as.ddmatrix(x)
  
  test(paste0("Random test #", i), {
    a <- na.exclude(x)
    b <- as.matrix(na.exclude(dx))
  }, check.attributes=FALSE)
}

for (i in 1:3){
  x <- matrix(na_maker(100), 100)
  dx <- as.ddmatrix(x, bldim=100)
  
  test(paste0("Random test (big bldim) #", i), {
    a <- na.exclude(x)
    b <- as.matrix(na.exclude(dx))
  }, check.attributes=FALSE)
}

collect()


finalize()

