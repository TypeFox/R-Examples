## ----, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = FALSE, comment = "#", message=FALSE)

## ----, echo=FALSE--------------------------------------------------------
require(Matrix)
require(trustOptim)
N <- 6
k <- 2
nv1 <- (N+1)*k
nels1 <- nv1^2
nnz1 <- (N+1)*k^2 + 2*N*k^2
nnz1LT <- (N+1)*k*(k+1)/2 + N*k^2
Q <- 1000
nv2 <- (Q+1)*k
nels2 <- nv2^2
nnz2 <- (Q+1)*k^2 + 2*Q*k^2
nnz2LT <- (Q+1)*k*(k+1)/2 + Q*k^2
options(scipen=999)

## ----, echo=FALSE--------------------------------------------------------
M <- as(kronecker(diag(N),matrix(1,k,k)),"lMatrix")
M <- rBind(M, Matrix(TRUE,k,N*k))
M <- cBind(M, Matrix(TRUE, k*(N+1), k))
print(M)

## ------------------------------------------------------------------------
set.seed(123)
data(binary)
str(binary)
N <- length(binary$Y)
k <- NROW(binary$X)
nvars <- as.integer(N*k + k)
start <- rnorm(nvars) ## random starting values
priors <- list(inv.Sigma = rWishart(1,k+5,diag(k))[,,1],
               inv.Omega = diag(k))

## ------------------------------------------------------------------------

opt <- trust.optim(start, fn=binary.f,
                   gr = binary.grad,  
                   hs = binary.hess,
                   method = "Sparse",
                   control = list(
                       start.trust.radius=5,
                       stop.trust.radius = 1e-7,
                       prec=1e-7,
                       report.precision=1L,
                       maxit=500L,
                       preconditioner=1L,
                       function.scale.factor=-1
                   ),
                   data=binary, priors=priors
                   )

