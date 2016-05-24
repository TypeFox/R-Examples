## ----setup1, echo=FALSE, cache=FALSE----------------------------------------------------
suppressPackageStartupMessages(require(Matrix))
knitr::render_sweave()
knitr::opts_chunk$set(prompt=TRUE, cache=TRUE)
options(replace.assign=TRUE, width=90,prompt="R> ")

## ----setup2, echo=FALSE-----------------------------------------------------------------
N <- 5
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

## ----blockarrow, echo=FALSE-------------------------------------------------------------
Mat <- as(kronecker(diag(N),matrix(1,k,k)),"sparseMatrix")
Mat <- rBind(Mat, Matrix(1,k,N*k))
Mat <- cBind(Mat, Matrix(1, k*(N+1), k))
printSpMatrix(as(Mat,"nMatrix"))

## ----banded, echo=FALSE-----------------------------------------------------------------
Mat <- kronecker(Matrix(1,k,k), diag(N))
Mat <- rBind(Mat, Matrix(1,k,N*k))
Mat <- cBind(Mat, Matrix(1, k*(N+1), k))
printSpMatrix(as(Mat,"nMatrix"))

## ---------------------------------------------------------------------------------------
library("sparseHessianFD")
bd <- kronecker(diag(3), matrix(TRUE,2,2))
Mat <- as(bd, "nMatrix")
printSpMatrix(tril(Mat))
mc <- Matrix.to.Coord(tril(Mat))
mc

## ---------------------------------------------------------------------------------------
pattern <- sparseMatrix(i=mc$rows, j=mc$cols)
printSpMatrix(pattern)

## ----eval=FALSE-------------------------------------------------------------------------
#  obj <- sparseHessianFD(x, fn, gr, rows, cols, ...)

## ----eval=FALSE-------------------------------------------------------------------------
#  f <- obj$fn(x)          ## returns numeric
#  df <- obj$gr(x)         ## returns numeric vector
#  hess <- obj$hessian(x)  ## returns dgCMatrix
#  fngr <- obj$fngr(x)     ## returns list
#  fngrhs <- obj$fngrhs(x) ## returns list

## ----binaryInit-------------------------------------------------------------------------
set.seed(123)
data("binary")
str(binary)
N <- length(binary[["Y"]])
k <- NROW(binary[["X"]])
T <- binary[["T"]]
nvars <- as.integer(N*k + k)
priors <- list(inv.Sigma = rWishart(1,k+5,diag(k))[,,1],
               inv.Omega = diag(k))

## ----trueValues-------------------------------------------------------------------------
P <- rnorm(nvars)
order.row <- FALSE
true.f <- binary.f(P, binary, priors, order.row=order.row)
true.grad <- binary.grad(P, binary, priors, order.row=order.row)
true.hess <- binary.hess(P, binary, priors, order.row=order.row)

## ----binaryRowsCols---------------------------------------------------------------------
pattern <- Matrix.to.Coord(tril(true.hess))
str(pattern)

## ----usingSparseHessianFD---------------------------------------------------------------
obj <- sparseHessianFD(P, fn=binary.f, gr=binary.grad,
       rows=pattern[["rows"]], cols=pattern[["cols"]],
       data=binary, priors=priors, order.row=order.row)
f <- obj$fn(P)
all.equal(f, true.f)
gr <- obj$gr(P)
all.equal(gr, true.grad)
hs <- obj$hessian(P)
all.equal(hs, true.hess)

