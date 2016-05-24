library(pbdTEST)
settings(mpi=TRUE)

.BLDIM <- 2
comm.set.seed(seed=1234, diff=FALSE)


### --------------------------------------
module("Reductions")

M <- 250
N <- 250
# M<- N<- 10
BL <- 4

x <- matrix(rnorm(M*N, 10, 100), M, N)
dx <- as.ddmatrix(x, BL)
x[1,1] <- x[1,1] + 1
dy <- as.ddmatrix(x, BL)

submodule("any, all, ==")

test("all(dx==dx)", {
  a <- all(dx==dx)
  b <- TRUE
})

test("any(dx==dx)", {
  a <- any(dx==dx)
  b <- TRUE
})

test("!all(dx==dB)", {
  a <- !all(dx==dy)
  b <- TRUE
})

test("any(dx==dB)", {
  a <- any(dx==dy)
  b <- TRUE
})


collect()





# ---------------------------------------------------
# Tests
# ---------------------------------------------------

tests <- function(.)
{
  rs1 <- rowSums(A)
  rs2 <- as.vector(rowSums(dx))
  comm.print(all.equal(rs2, rs2))

  out1 <- colSums(A)
  out2 <- as.vector(colSums(dx))
  comm.print(all.equal(out1, out2))
  
  rs1 <- rowMeans(A)
  rs2 <- as.vector(rowMeans(dx))
  comm.print(all.equal(rs1, rs2))
#  if (!is.logical(all.equal(rs1, rs2))){
#    comm.print(rs1)
#    comm.print(rs2)
#  }
  
  out1 <- colMeans(A)
  out2 <- as.vector(colMeans(dx))
  comm.print(all.equal(out1, out2))
#  if (!is.logical(all.equal(out1, out2))){
#    comm.print(out1)
#    comm.print(out2)
#  }
  
  out1 <- sum(A)
  out2 <- sum(dx)
  comm.print(all.equal(out1, out2))
  
  out1 <- prod(A)
  out2 <- prod(dx)
  comm.print(all.equal(out1, out2))

  out1 <- diag(A)
  out2 <- diag(dx)
  comm.print(all.equal(out1, out2))
#  if (!is.logical(all.equal(out1, out2))){
#    comm.print(out1)
#    comm.print(out2)
#  }
  
  out1 <- mean(A)
  out2 <- mean(dx)
  comm.print(all.equal(out1, out2))
  
  out1 <- apply(A, MARGIN=2, FUN=sd)
  out2 <- as.vector(sd(dx))
  comm.print(all.equal(out1, out2))
}



comm.print("-------Reductions-------")
comm.print("       Square")

A <- matrix(rnorm(M*N, 10, 100), M, N)
dx <- as.ddmatrix(A, BL)
tests()

comm.print("       Column")

A <- matrix(rnorm(M*1, 10, 100), M, 1)
dx <- as.ddmatrix(A, BL)
tests()

comm.print("       Row")

A <- matrix(rnorm(1*N, 10, 100), 1, N)
dx <- as.ddmatrix(A, BL)
tests()

finalize()
