## ----, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#", message=FALSE)
options(digits=4, scipen=0)			       

## ------------------------------------------------------------------------
require(trustOptim)
require(Matrix)
f <- function(V) {

    N <- length(V)/2
    x <- V[seq(1,2*N-1,by=2)]
    y <- V[seq(2,2*N,by=2)]
    return(sum(100*(x^2-y)^2+(x-1)^2))
}

df <- function(V) {
    N <- length(V)/2
    x <- V[seq(1,2*N-1,by=2)]
    y <- V[seq(2,2*N,by=2)]

    t <- x^2-y
    dxi <- 400*t*x+2*(x-1)
    dyi <- -200*t
    return(as.vector(rbind(dxi,dyi)))
 }

hess <- function(V) {

    N <- length(V)/2
    x <- V[seq(1,2*N-1,by=2)]
    y <- V[seq(2,2*N,by=2)]
    d0 <- rep(200,N*2)
    d0[seq(1,(2*N-1),by=2)] <- 1200*x^2-400*y+2
    d1 <- rep(0,2*N-1)
    d1[seq(1,(2*N-1),by=2)] <- -400*x

    H <- bandSparse(2*N,
                    k=c(-1,0,1),
                    diagonals=list(d1,d0,d1),
                    symmetric=FALSE,
                    giveCsparse=TRUE)
    return(drop0(H))
}

## ------------------------------------------------------------------------
set.seed(1234)
N <- 3
start <- as.vector(rnorm(2*N, -1, 3))

## ------------------------------------------------------------------------
opt <- trust.optim(start, fn=f, gr=df, hs=hess, method="Sparse")

## ------------------------------------------------------------------------
opt

## ------------------------------------------------------------------------
opt1 <- trust.optim(start, fn=f, gr=df, hs=hess, method="Sparse",
      control=list(preconditioner=1, report.freq=10))

## ------------------------------------------------------------------------
opt.bfgs <- trust.optim(start, fn=f, gr=df, method="BFGS", control=list(maxit=5000, report.freq=10))
opt.bfgs

## ------------------------------------------------------------------------
opt.sr1 <- trust.optim(start, fn=f, gr=df, method="SR1", control=list(maxit=5000, report.freq=10))
opt.sr1

