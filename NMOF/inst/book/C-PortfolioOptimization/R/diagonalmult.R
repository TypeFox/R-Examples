## diagonalmult.R -- version 2010-11-27

# set size of matrix
N <- 1000L

# create matrices
D <- diag(runif(N))
C <- array(runif(N * N), dim = c(N, N))

# compute product / compare time
system.time(Z1 <- D %*% C %*% D)
system.time(Z2 <- outer(diag(D), diag(D)) * C)

# check difference between matrices
max(abs(Z1 - Z2)) # ... or use all.equal(Z1,Z2)

# ... or with the Matrix package
require(Matrix)
D2 <- Diagonal(x = diag(D))
system.time(Z3 <- D2 %*% C %*% D2)
system.time(Z4 <- outer(diag(D2), diag(D2)) * C)
stopifnot(all.equal(as.numeric(Z1),as.numeric(Z4)))