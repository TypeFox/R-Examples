# ################################################
# ------------------------------------------------
# ScaLAPACK
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
# ScaLAPACK
# ---------------------------------------------------

comm.print("-------ScaLAPACK solve()-------")

A <- matrix(rnorm(M*N, 10, 100), M, N)
dA <- as.ddmatrix(A, BL)

B <- matrix(rnorm(M*3, 10, 100), M, 3)
dB <- as.ddmatrix(B, BL)

out1 <- solve(A, B)
out2 <- as.matrix(solve(dA, dB))
comm.print(all.equal(out1, out2))

out1 <- solve(A)
out2 <- as.matrix(solve(dA))
comm.print(all.equal(out1, out2))



comm.print("-------ScaLAPACK svd()-------")
comm.print("       Square")

A <- matrix(rnorm(M*N, 10, 100), M, N)
dA <- as.ddmatrix(A, BL)

  out1 <- La.svd(A)$d
  out2 <- La.svd(dA)$d
  comm.print(all.equal(out1, out2))

  out1 <- La.svd(A, nu=0)$d
  out2 <- La.svd(dA, nu=0)$d
  comm.print(all.equal(out1, out2))

  out1 <- La.svd(A, nv=0)$d
  out2 <- La.svd(dA, nv=0)$d
  comm.print(all.equal(out1, out2))

  out1 <- La.svd(A, nv=1)$d
  out2 <- La.svd(dA, nv=1)$d
  comm.print(all.equal(out1, out2))

  out1 <- svd(A)$d
  out2 <- svd(dA)$d
  comm.print(all.equal(out1, out2))

  out1 <- svd(A, nu=0)$d
  out2 <- svd(dA, nu=0)$d
  comm.print(all.equal(out1, out2))

  out1 <- svd(A, nv=0)$d
  out2 <- svd(dA, nv=0)$d
  comm.print(all.equal(out1, out2))



comm.print("       Column")

A <- matrix(rnorm(M*1, 10, 100), M, 1)
dA <- as.ddmatrix(A, BL)

out1 <- svd(A)$d
out2 <- svd(dA)$d
comm.print(all.equal(out1, out2))

out1 <- svd(A, nu=0)$d
out2 <- svd(dA, nu=0)$d
comm.print(all.equal(out1, out2))

out1 <- svd(A, nv=0)$d
out2 <- svd(dA, nv=0)$d
comm.print(all.equal(out1, out2))



comm.print("       Row")

A <- matrix(rnorm(1*N, 10, 100), 1, N)
dA <- as.ddmatrix(A, BL)

out1 <- svd(A)$d
out2 <- svd(dA)$d
comm.print(all.equal(out1, out2))

out1 <- svd(A, nu=0)$d
out2 <- svd(dA, nu=0)$d
comm.print(all.equal(out1, out2))

out1 <- svd(A, nv=0)$d
out2 <- svd(dA, nv=0)$d
comm.print(all.equal(out1, out2))



comm.print("-------ScaLAPACK chol()-------")

A <- crossprod(matrix(rnorm(N*N, 10, 100), N, N))
dA <- as.ddmatrix(A, BL)

cholA <- chol(A)
choldA <- chol(dA)

out1 <- chol(A)
out2 <- as.matrix(chol(dA))
comm.print(all.equal(out1, out2))

out1 <- chol2inv(x=cholA)
out2 <- as.matrix(chol2inv(choldA))
comm.print(all.equal(out1, out2))



comm.print("-------ScaLAPACK lu()-------")

N <- 25
A <- matrix(rnorm(N*N, 10, 100), N, N)
dA <- as.ddmatrix(A, BL)

out2 <- as.matrix( lu(dA) )

suppressPackageStartupMessages(library(Matrix))
out1 <- matrix(lu(A)@x, nrow=N, ncol=N)
comm.print(all.equal(out1, out2))


finalize()
