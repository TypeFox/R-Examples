#' Quadratic monotone spline basis function for given knots.
#' 
#' Calculate basis functions for monotone quadratic splines.
#' 
#' @param xvec Vector at which to evaluate the basis functions.
#' @param tvec Vector of spline knots: lower endpoint, interior knot, upper endpoint.
#' @param intercept Logical; should an intercept be included or not?
#' @keywords multivariate
#' @export ispline
ispline <- function (xvec, tvec = c(0, 0.5, 1), intercept = TRUE) {
  stopifnot(length(tvec) == 3)
    n <- length(xvec)
    i1 <- function(x) {
        n <- length(x)
        out <- rep(0, n)
        for (i in 1:n) {
            if (x[i] >= tvec[1] & x[i] < tvec[2]) 
                out[i] <- (2 * tvec[2] * (x[i] - tvec[1]) - (x[i]^2 - 
                  tvec[1]^2))/(tvec[2] - tvec[1])^2
            if (x[i] < tvec[1]) 
                out[i] <- 0
            if (x[i] >= tvec[2]) 
                out[i] <- 1
        }
        out
    }
    i2 <- function(x) {
        n <- length(x)
        out <- rep(0, n)
        for (i in 1:n) {
            if (x[i] >= tvec[1] & x[i] < tvec[2]) 
                out[i] <- (x[i] - tvec[1])^2/((tvec[2] - tvec[1]) * 
                  (tvec[3] - tvec[1]))
            if (x[i] >= tvec[2] & x[i] < tvec[3]) 
                out[i] <- (tvec[2] - tvec[1])/(tvec[3] - tvec[1]) + 
                  (2 * tvec[3] * (x[i] - tvec[2]) - (x[i]^2 - 
                    tvec[2]^2))/((tvec[3] - tvec[2]) * (tvec[3] - 
                    tvec[1]))
            if (x[i] < tvec[1]) 
                out[i] <- 0
            if (x[i] >= tvec[3]) 
                out[i] <- 1
        }
        out
    }
    i3 <- function(x) {
        n <- length(x)
        out <- rep(0, n)
        for (i in 1:n) {
            if (x[i] >= tvec[2] & x[i] < tvec[3]) 
                out[i] <- (x[i] - tvec[2])^2/(tvec[3] - tvec[2])^2
            if (x[i] < tvec[2]) 
                out[i] <- 0
            if (x[i] >= tvec[3]) 
                out[i] <- 1
        }
        out
    }
    out <- matrix(0, ncol = 3, nrow = n)
    out[, 1] <- i1(xvec)
    out[, 2] <- i2(xvec)
    out[, 3] <- i3(xvec)
    colnames(out) <- paste0("m", 1:3)
    if (intercept) {
        out <- cbind(1, out)
        colnames(out)[1] <- "(Intercept)"
    }
    out
}