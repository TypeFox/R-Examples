## ----, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = FALSE, comment = "#", message=FALSE)
options(digits=4)

## ----, echo=FALSE--------------------------------------------------------
require(sparseMVN)
require(Matrix)
require(mvtnorm)
N <- 5
k <- 2
p <- 2
nv1 <- N*k+p
nels1 <- nv1^2
nnz1 <- N*k^2 + 2*p*N*k + p^2
nnz1LT <- N*k*(k+1)/2  + p*N*k + p*(p+1)/2
Q <- 1000
nv2 <- Q*k+p
nels2 <- nv2^2
nnz2 <- Q*k^2 + 2*p*Q*k + p^2
nnz2LT <- Q*k*(k+1)/2 + p*Q*k + p*(p+1)/2
options(scipen=999)

## ----, echo=FALSE--------------------------------------------------------
M <- as(kronecker(diag(N),matrix(1,k,k)),"lMatrix")
M <- rBind(M, Matrix(TRUE,p,N*k))
M <- cBind(M, Matrix(TRUE, k*N+p, p))
print(M)

## ----, eval=FALSE--------------------------------------------------------
#  rmvn.sparse(ndraws, mu, CH, prec=TRUE)
#  dmvn.sparse(x, mu, CH, prec=TRUE)

## ----, collapse=TRUE-----------------------------------------------------
require(Matrix)
set.seed(123)
N <- 5  ## number of blocks in sparse covariance matrix
p <- 2 ## size of each block
k <- 2  ##
R <- 10
    
## mean vector
mu <- seq(-3,3,length=k*N+p)

## build random block-arrow covariance/precision matrix for test
Q1 <- tril(kronecker(diag(N), Matrix(seq(0.1,1.1,length=k*k),k,k)))
Q2 <- Matrix(rnorm(N*k*p), p, N*k)
Q3 <- tril(0.2*diag(p))
Sigma <- rBind(cBind(Q1, Matrix(0, N*k, p)), cBind(Q2, Q3))
Sigma <- Matrix::tcrossprod(Sigma)
class(Sigma)

## ----, collapse=TRUE-----------------------------------------------------
chol.Sigma <- Matrix::Cholesky(Sigma)  ## creates a dCHMsimpl object
x <- rmvn.sparse(R, mu, chol.Sigma, prec=FALSE)

## ----, collapse=TRUE-----------------------------------------------------
d <- dmvn.sparse(x, mu, chol.Sigma, prec=FALSE)
str(d)

