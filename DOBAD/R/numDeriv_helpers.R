### adapting the functions from CRAN package numDeriv so that they work
### for differentiation from below.  Note: an enormous amount of accuracy is
### apparently lost because of the non-cancellation of odd terms.  Not sure
### if there is any sort of work-around.





##sides: +1 from above, -1 from below
## might be able to get better approx by allowing some to be 2-sided..
genDoneSided <- function (func, x, sides, method = "Richardson",
                          method.args = list(eps = 1e-04, d = 1e-04, zero.tol = sqrt(.Machine$double.eps/7e-07), r = 4, 
                            v = 2), ...) 
{
  if (method != "Richardson") 
    stop("method not implemented.")
  if (length(x) != sum(abs(sides)))
      stop("genDoneSided: length(x) != length(sides) or bad arg passed for sides.");
  args <- list(eps = 1e-04, d = 1e-04, zero.tol = sqrt(.Machine$double.eps/7e-07), 
               r = 4, v = 2)
  args[names(method.args)] <- method.args
  eps <- args$eps
  d <- args$d
  r <- args$r
  v <- args$v
  if (v != 2) 
    stop("The current code assumes v is 2 (the default).") #great...
  f0 <- func(x, ...)
  p <- length(x)
  ##h0 <- abs(d * x) + eps * (abs(x) < args$zero.tol) # standard code
  h0 <- abs(d * x)*(abs(x) >= args$zero.tol)  + eps*(abs(x) < args$zero.tol)
  D <- matrix(0, length(f0), (p * (p + 3))/2)
  Daprox <- matrix(0, length(f0), r)
  Hdiag <- matrix(0, length(f0), p)
  Haprox <- matrix(0, length(f0), r)
  for (i in 1:p) {
    h <- h0
    for (k in 1:r) {
      ##f1 <- func(x + (i == (1:p)) * h, ...)
      f2 <- func(x + (sides*(i == (1:p))) * h, ...)
      f1 <- f0;
      f3 <- func(x + (sides* (i == (1:p))) *2* h, ...)
      ##Daprox[, k] <- (f1 - f2)/(2 * h[i])
      Daprox[, k] <- (f2-f1)/(h[i]*sides[i])
      ##Haprox[, k] <- (f1 - 2 * f0 + f2)/h[i]^2
      Haprox[, k] <- (f1 - 2 * f2 + f3)/h[i]^2
      h <- h/v
      NULL
    }
    for (m in 1:(r - 1)) for (k in 1:(r - m)) {
      ##             Daprox[, k] <- (Daprox[, k + 1] * (4^m) - Daprox[, 
      ##                 k])/(4^m - 1)
      ##             Haprox[, k] <- (Haprox[, k + 1] * (4^m) - Haprox[, 
      ##                 k])/(4^m - 1)
      Daprox[, k] <- (Daprox[, k + 1] * (2^m) -
                      Daprox[,k])/(2^m - 1)
      Haprox[, k] <- (Haprox[, k + 1] * (2^m) -
                      Haprox[,k])/(2^m - 1)
      NULL
    }
    D[, i] <- Daprox[, 1]
    Hdiag[, i] <- Haprox[, 1]
    NULL
  }
  u <- p
  for (i in 1:p) {
    for (j in 1:i) {
      u <- u + 1
      if (i == j) {
        D[, u] <- Hdiag[, i]
        NULL
      }
      else {
        h <- h0
        for (k in 1:r) {
          f1  <- f0;
          f4  <- func(x + (sides*(i == (1:p))) * h + (sides*(j == (1:p)))*h, ...)
          f01 <- func(x + (sides*(i == (1:p))) * h , ...)
          f02 <- func(x + (sides*(j == (1:p))) * h , ...)
          ###f3 <- func(x + (sides*(i == (1:p))) * h + (sides* (j == (1:p)))*h, ...)
          Daprox[, k] <- (f1 - f01 -f02 + f4)/(sides[i]*sides[j]*h[i]*h[j])
        }
        for (m in 1:(r - 1)) for (k in 1:(r - m)) {
          Daprox[, k] <- (Daprox[, k + 1] * (2^m) -
                          Daprox[, k])/(2^m - 1)
          NULL
        }
        D[, u] <- Daprox[, 1]
        NULL
      }
    }
  }
  D <- list(D = D, p = length(x), f0 = f0, func = func, x = x, 
            d = d, method = method, method.args = args)
  class(D) <- "Darray"
  invisible(D)
}


hessianOneSided <- function (func, x, sides, method = "Richardson", method.args = list(eps = 1e-04, 
                                                            d = 1e-4, zero.tol = sqrt(.Machine$double.eps/7e-07), r = 4, 
                                                            v = 2), ...) 
{
  if (method != "Richardson") 
    stop("method not implemented.")
  if (1 != length(func(x, ...))) 
    stop("hessian.default assumes a scalar real valued function.")
  D <- genDoneSided(func, x,sides=sides, method = method, method.args = method.args, 
                 ...)$D
  if (1 != nrow(D)) 
    stop("BUG! should not get here.")
  H <- diag(NA, length(x))
  u <- length(x)
  for (i in 1:length(x)) {
    for (j in 1:i) {
      u <- u + 1
      H[i, j] <- D[, u]
      H[j, i] <- D[, u]
    }
  }
  ##H <- H + t(H)
  ##diag(H) <- diag(H)/2
  H
}

## ###### i don't like how they use epsilon and d, but ...
## ## ###this is from the numDeriv package. however,
## genD.2 <- function (func, x, method = "Richardson",
##                   method.args=list(eps = 1e-04, 
##                   d = 1e-04, zero.tol = sqrt(.Machine$double.eps/7e-07),
##                   r = 4, v = 2),
##                   ...) 
## {
##   if (method != "Richardson") 
##     stop("method not implemented.")
##   ###this is the original code ...
##   ##i dont understand the point
##   args <- list(eps = 1e-04, d = 1e-04,
##   zero.tol = sqrt(.Machine$double.eps/7e-07), 
##                r = 4, v = 2)
##   args[names(method.args)] <- method.args
##   print( paste("method.args", method.args))
##   print("args"); print(args);
##   eps <- args$eps
##   d <- args$d
##   r <- args$r
##   v <- args$v
##   if (v != 2) {
##     print(paste("v is " , v));
##         stop("The current code assumes v is 2 (the default).")
##       }
##     f0 <- func(x, ...)
##     p <- length(x)
##   h0 <- abs(d * x) + eps * (abs(x) < args$zero.tol)
##   D <- matrix(0, length(f0), (p * (p + 3))/2)
##     Daprox <- matrix(0, length(f0), r)
##     Hdiag <- matrix(0, length(f0), p)
##     Haprox <- matrix(0, length(f0), r)
##     for (i in 1:p) {
##         h <- h0
##         for (k in 1:r) {
##             f1 <- func(x + (i == (1:p)) * h, ...)
##             f2 <- func(x - (i == (1:p)) * h, ...)
##             Daprox[, k] <- (f1 - f2)/(2 * h[i])
##             Haprox[, k] <- (f1 - 2 * f0 + f2)/h[i]^2
##             h <- h/v
##             NULL
##         }
##         for (m in 1:(r - 1)) for (k in 1:(r - m)) {
##             Daprox[, k] <- (Daprox[, k + 1] * (4^m) - Daprox[, 
##                 k])/(4^m - 1)
##             Haprox[, k] <- (Haprox[, k + 1] * (4^m) - Haprox[, 
##                 k])/(4^m - 1)
##             NULL
##         }
##         D[, i] <- Daprox[, 1]
##         Hdiag[, i] <- Haprox[, 1]
##         NULL
##     }
##     u <- p
##     for (i in 1:p) {
##         for (j in 1:i) {
##             u <- u + 1
##             if (i == j) {
##                 D[, u] <- Hdiag[, i]
##                 NULL
##             }
##             else {
##                 h <- h0
##                 for (k in 1:r) {
##                   f1 <- func(x + (i == (1:p)) * h + (j == (1:p)) * 
##                     h, ...)
##                   f2 <- func(x - (i == (1:p)) * h - (j == (1:p)) * 
##                     h, ...)
##                   Daprox[, k] <- (f1 - 2 * f0 + f2 - Hdiag[, 
##                     i] * h[i]^2 - Hdiag[, j] * h[j]^2)/(2 * h[i] * 
##                     h[j])
##                   h <- h/v
##                 }
##                 for (m in 1:(r - 1)) for (k in 1:(r - m)) {
##                   Daprox[, k] <- (Daprox[, k + 1] * (4^m) - Daprox[, 
##                     k])/(4^m - 1)
##                   NULL
##                 }
##                 D[, u] <- Daprox[, 1]
##                 NULL
##             }
##         }
##     }
##     D <- list(D = D, p = length(x), f0 = f0, func = func, x = x, 
##         d = d, method = method, method.args = args)
##     class(D) <- "Darray"
##     invisible(D)
## }


## ###no change from numDeriv::hessian but
## ### now may have complications with genD if i load numDeriv.
## hessian <- function (func, x, method = "Richardson", method.args = list(eps = 1e-04, 
##     d = 0.1, zero.tol = sqrt(.Machine$double.eps/7e-07), r = 4, 
##     v = 2), ...) 
## {
##     if (method != "Richardson") 
##         stop("method not implemented.")
##     if (1 != length(func(x, ...))) 
##         stop("hessian.default assumes a scalar real valued function.")
##     D <- genD(func, x, method = method, method.args = method.args, 
##         ...)$D
##     if (1 != nrow(D)) 
##         stop("BUG! should not get here.")
##     H <- diag(NA, length(x))
##     u <- length(x)
##     for (i in 1:length(x)) {
##         for (j in 1:i) {
##             u <- u + 1
##             H[i, j] <- D[, u]
##         }
##     }
##     H <- H + t(H)
##     diag(H) <- diag(H)/2
##     H
## }











#########################
###########more dummy versions of differntiation functions:
##########################

## Differentiate either a real or complex function, where a possible complex
## argument is a function of one real variable

num.deriv = function(ftn, var, delta=0.001,...){
  return((ftn(var+delta,...) - ftn(var-delta,...))/(2*delta))
}


#### not used; prefer modified functions from 'numDeriv' package.
num.deriv2 <- function(ftn, var1,var2, delta=0.00001){
  (ftn(c(var1+delta, var2+delta)) + ftn(c(var1-delta, var2-delta)) -
    ftn(c(var1+delta, var2-delta)) - ftn(c(var1-delta, var2+delta)))/(4*delta^2);
}

#2nd deriv of f:R -> R
num.2deriv <- function(ftn, var, delta=0.00001,...){
  num.deriv(ftn=function(x){num.deriv(ftn=ftn, var=x, delta=delta,...);},
            var=var,
            delta=delta);
}
