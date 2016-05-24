rlogcon <- function(n, x0){

## ===================================================
## Simulate from $\hat f_n$ and $\hat f_n^*$ as described
## in Duembgen and Rufibach (2009, Section 3)
## ===================================================

## Set parameters
res <- logConDens(x0, xgrid = NULL, smoothed = FALSE, print = FALSE, gam = NULL, xs = NULL)
phi <- res$phi
f.smoothed <- evaluateLogConDens(res$x, res, which = 4, gam = NULL)[, "smooth.density"] 

## compute gam
VarFn <- LocalVariance(x = res$x, w = res$w, phi = phi)
gam <- sqrt(res$sig ^ 2 - VarFn)

## Simulate from $\hat f_n$ and $\hat f_n^*$
X <- rep(NA, length = n)
U <- runif(n)
Z <- rnorm(n)
X <- quantilesLogConDens(U, res)[, "quantile"]

X_star <- X + gam * Z

res <- list("X" = X, "X_star" = X_star, "U" = U, "Z" = Z, "f" = exp(phi), "f.smoothed" = f.smoothed, "x" = res$x, "w" = res$w)
return(res)
}




