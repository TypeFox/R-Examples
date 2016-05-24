# ################################################
# ------------------------------------------------
# Arithmetic
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
#M <- N <- 10
BL <- 4

comm.set.seed(seed=1234, diff=F)


tol <- 1e-8

# ---------------------------------------------------
# Tests
# ---------------------------------------------------

tests <- function(.)
{
  out1 <- 3*A-A + A
  out2 <- as.matrix(3*dA-dA+dA)
  comm.print(all.equal(out1, out2))
  
  out1 <- sum(A+1)
  out2 <- sum(dA+1)
  comm.print(all.equal(out1, out2))

  out1 <- 4+A+1 + A/10
  out2 <- as.matrix(  4+dA+1 + dA/10  )
  comm.print(all.equal(out1, out2))

  out1 <- prod(A/100)
  out2 <- prod(dA/100)
  comm.print(all.equal(out1, out2))
  
  out1 <- log(sqrt(abs(A)))
  out2 <- as.matrix(log(sqrt(abs(dA))))
  comm.print(all.equal(out1, out2))
}

# ---------------------------------------------------
# Arithmetic
# ---------------------------------------------------


M <- N <- 10

comm.print("-------Arithmetic-------")
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

comm.print("-------Extra-------")

A <- matrix(rnorm(M*N, 10, 100), M, N)
dA <- as.ddmatrix(A, BL)
vec <- 1:5

comm.print("       +")
  # +
  out1 <- A+vec
  out2 <- as.matrix(dA+vec)
  comm.print(all( abs(out1-out2) < tol))

  out1 <- vec+A
  out2 <- as.matrix(vec+dA)
  comm.print(all( abs(out1-out2) < tol))

comm.print("       -")
  # - 
  out1 <- A-vec
  out2 <- as.matrix(dA-vec)
  comm.print(all( abs(out1-out2) < tol))

  out1 <- vec-A
  out2 <- as.matrix(vec-dA)
  comm.print(all( abs(out1-out2) < tol))

comm.print("       *")
  # * 
  out1 <- A*vec
  out2 <- as.matrix(dA*vec)
  comm.print(all( abs(out1-out2) < tol))

  out1 <- vec*A
  out2 <- as.matrix(vec*dA)
  comm.print(all( abs(out1-out2) < tol))

comm.print("       /")
  # /
  out1 <- A/vec
  out2 <- as.matrix(dA/vec)
  comm.print(all( abs(out1-out2) < tol))

  out1 <- vec/A
  out2 <- as.matrix(vec/dA)
  comm.print(all( abs(out1-out2) < tol))

comm.print("       ^")
  # ^
  out1 <- A^vec
  out2 <- as.matrix(dA^vec)
  comm.print(all( abs(out1-out2) < tol))

comm.print("       %%")
  # %%
  vec <- c(-3, -2, -1, 1, 2)

  o1 <- A %% vec
  o2 <- as.matrix( dA %% vec )
  comm.print(all( abs(o1-o2) <tol))

  o1 <- vec %% A
  o2 <- as.matrix( vec %% dA )
  comm.print(all( abs(o1-o2) <tol))

  vec <- 0

  o1 <- A %% vec
  o2 <- as.matrix( dA %% vec )
  comm.print(all(is.na(o1)) && all(is.na(o2)))

  o1 <- vec %% A
  o2 <- as.matrix( vec %% dA )
  comm.print(all( abs(o1-o2) <tol))

comm.print("       %/%")
  # /
  vec <- 1:5
  
  out1 <- A%/%vec
  out2 <- as.matrix(dA%/%vec)
  comm.print(all( abs(out1-out2) < tol))

  out1 <- vec%/%A
  out2 <- as.matrix(vec%/%dA)
  comm.print(all( abs(out1-out2) < tol))


finalize()
