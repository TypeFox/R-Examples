# ################################################
# ------------------------------------------------
# PBLAS
# ------------------------------------------------
# ################################################

# For each test, script returns TRUE if the test was successful 
# (produced the correct value), and returns FALSE if the test was
# unsuccessful (produced the incorrect value).


suppressPackageStartupMessages(library(pbdDMAT, quiet=T))

init.grid()

comm.print <- function(x) pbdMPI::comm.print(x, quiet=T)

M <- 250
N <- 250
BL <- 4

comm.set.seed(seed=1234, diff=F)

tol <- 1e-8

# ---------------------------------------------------
# Tests
# ---------------------------------------------------

tests <- function(.)
{
#  out1 <- t(A) %*% A
#  out2 <- as.matrix(t(dA) %*% dA)
#  comm.print(all.equal(out1, out2))
  
  out1 <- crossprod(A)
  out2 <- as.matrix(crossprod(dA))
  comm.print(all.equal(out1, out2))
  
  out1 <- tcrossprod(A)
  out2 <- as.matrix(tcrossprod(dA))
  comm.print(all.equal(out1, out2))
  
  out1 <- crossprod(A, A)
  out2 <- as.matrix(crossprod(dA, dA))
  comm.print(all.equal(out1, out2))
  
}

# ---------------------------------------------------
# PBLAS
# ---------------------------------------------------

comm.print("-------PBLAS-------")
comm.print("       Square")

A <- matrix(rnorm(M*N, 10, 100), M, N)
dA <- as.ddmatrix(A, BL)
tests()

comm.print("       Column")

A <- matrix(rnorm(M*1, 10, 100), M, 1)
dA <- as.ddmatrix(A, BL)
tests()

comm.print("       Row")

A <- matrix(rnorm(1*N, 10, 100), 1, N)
dA <- as.ddmatrix(A, BL)
tests()

finalize()
