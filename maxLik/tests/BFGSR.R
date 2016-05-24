### BFGSR-related tests

## 1. Test maximization algorithm for convex regions
## 
## Optimize quadratic form t(D) %*% W %*% D with p.d. weight matrix
## (ie unbounded problems).
## All solutions should go to large values with a message about successful convergence
set.seed(0)
options(digits=4)
quadForm <- function(D) {
   return( - t(D - (1:N) ) %*% W %*% ( D - (1:N) ) )
}
N <- 3

# round gradients to increase reproducibility of the accuracy
roundGradients <- function( object ) {
   object$gradient <- round( object$gradient, 3 )
   return( object )   
}
                           # 3-dimensional case
## a) test quadratic function t(D) %*% D
library(maxLik)
W <- diag(N)
D <- rep(1/N, N)
res <- maxBFGSR(quadForm, start=D)
res <- roundGradients( res )
summary(res)

## b) add noice to
W <- diag(N) + matrix(runif(N*N), N, N)
                           # diagonal weight matrix with some noise
D <- rep(1/N, N)
res <- maxBFGSR(quadForm, start=D, tol = 1e-10 )
res <- roundGradients( res )
summary(res)

## Next, optimize hat function in non-concave region.  Does not work well.
hat <- function(param) {
   ## Hat function.  Hessian negative definite if sqrt(x^2 + y^2) < 0.5
   x <- param[1]
   y <- param[2]
   exp(-(x-2)^2 - (y-2)^2)
}

hatNC <- maxBFGSR(hat, start=c(1,1), tol=0, reltol=0)
hatNC <- roundGradients( hatNC )
summary( hatNC )

# should converge to c(0,0).

## Test BFGSR with fixed parameters and equality constraints
## Optimize 3D hat with one parameter fixed (== 2D hat).
## Add an equality constraint on that
hat3 <- function(param) {
   ## Hat function.  Hessian negative definite if sqrt((x-2)^2 + (y-2)^2) < 0.5
   x <- param[1]
   y <- param[2]
   z <- param[3]
   exp(-(x-2)^2-(y-2)^2-(z-2)^2)
}
sv <- c(x=1,y=1,z=1)
## constraints: x + y + z = 8
A <- matrix(c(1,1,1), 1, 3)
B <- -8
constraints <- list(eqA=A, eqB=B)
hat3CF <- maxBFGSR(hat3, start=sv, constraints=constraints, fixed=3)
hat3CF <- roundGradients( hat3CF )
summary( hat3CF )
