##
##  n e l d e r _ m e a d . R  Nelder-Mead Minimization
##


nelder_mead <- function(x0, f, lb = NULL, ub = NULL, tol = 1e-10,
                        maxfeval = 20000, step = rep(1.0, length(x0)), ...) {
    if (!is.numeric(x0))
        stop("Argument 'x0' must be a numeric vector.",
              call. = FALSE)
    n <- length(x0)
    if (n == 1)
        stop("Do not use Nelder-Mead for univariate functions.",
              call. = FALSE)
    if (length(step) != n)
        stop("Argument 'step' must be of the same length as 'x0'.",
             call. = FALSE)

    if (!is.function(f) || is.null(f))
        stop("Function argument 'f' must be a valid R function.",
              call. = FALSE)
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    if (is.null(lb) && is.null(ub)) {
        result <- .nelmin(x0, f, tol, maxfeval, step)
    } else {
        result <- .nelminb(x0, f, lb, ub, tol, maxfeval, step)
    }

    return(result)
}


.nelmin <- function(x0, fn, tol, maxfeval, step) {
    n <- length(x0)
                        # Inputs:
    start  <- x0        # starting point
    reqmin <- tol       # terminating limit for the variance of function values
    # step <- step      # size and shape of the initial simplex
    kcount <- maxfeval  # maximum number of function evaluations.
    konvge <- kcount/100# convergence check is carried out every
                        # KONVGE iterations, >= 1

                        # Outputs:
    xmin   <- NA        # estimated minimum of the function
    ynewlo <- NA        # minimum value of the function
    icount <- 0         # number of function evaluations
    numres <- 0         # number of restarts, must be > 1
    ifault <- 0         # error indicator, 0, 1, 2

    # Constants for Nelder-Mead
    ccoeff <- 0.5
    ecoeff <- 2.0
    eps    <- 0.001
    rcoeff <- 1.0
    
    jcount <- konvge
    dn <- n
    nn <- n + 1
    dnn <- nn
    del <- 1.0
    rq <- reqmin * dn

    pbar   <- numeric(n)  # centroid
    pstar  <- numeric(n)
    p2star <- numeric(n)

    y <- numeric(n+1)
    p <- matrix(0, n, n+1)

    #-- Start the main loop ------------
    while ( TRUE ) {                    # outer while loop
        # Initial or restarted loop
        p[, nn] <- start
        y[nn] <- fn ( start )
        icount <- icount + 1

        for (j in 1:n) {
          x <- start[j]
          start[j] <- start[j] + step[j] * del
          p[, j] <- start
          y[j] <- fn ( start )
          icount <- icount + 1
          start[j] <- x
        } # simplex construction is complete

        # Find highest and lowest y values.
        ilo <- which.min(y)
        ylo <- y[ilo]

        iinner <- 0
        while ( icount < kcount ) {                # inner while loop
            iinner <- iinner + 1
            # indicate the vertex of the simplex to be replaced
            ihi <- which.max(y)
            ynewlo <- y[ihi]

            # Calculate the centroid of the simplex vertices
            #   excepting the vertex with Y value YNEWLO
            pbar <- rowSums(p[, -ihi]) / dn

            # Reflection through the centroid
            pstar <- pbar + rcoeff * ( pbar - p[,ihi] )
            ystar <- fn ( pstar )
            icount <- icount + 1

            # Successful reflection, so extension
            if ( ystar < ylo ) {
                p2star = pbar + ecoeff * ( pstar - pbar )
                y2star <- fn ( p2star )
                icount <- icount + 1

                # Check extension.
                if ( ystar < y2star ) {
                    p[, ihi] <- pstar
                    y[ihi] <- ystar

                # Retain extension or contraction.
                } else {
                    p[, ihi] <- p2star
                    y[ihi] <- y2star
                }

            #  No extension.
            } else {
                l <- sum(ystar < y)
                if ( l > 1 ) {
                    p[, ihi] <- pstar
                    y[ihi] <- ystar

                # Contraction on the Y(IHI) side of the centroid.
                } else if ( l == 0 ) {
                    p2star <- pbar + ccoeff * ( p[, ihi] - pbar )
                    y2star <- fn ( p2star )
                    icount <- icount + 1

                    #  Contract the whole simplex.
                    if ( y[ihi] < y2star ) {
                        for (j in 1:nn) {
                            p[, j] <- (p[, j] + p[, ilo]) / 2
                            xmin <- p[, j]
                            y[j] <- fn ( xmin )
                            icount <- icount + 1
                        }

                        ilo <- which.min(y)
                        ylo <- y[ilo]
                        next

                    # Retain contraction
                    } else {
                        p[, ihi] <- p2star
                        y[ihi] <- y2star
                    }

                #  Contraction on the reflection side of the centroid
                } else if ( l == 1 ) {
                    p2star <- pbar + ccoeff * ( pstar - pbar )
                    y2star <- fn ( p2star )
                    icount <- icount + 1

                    #  Retain reflection?
                    if ( y2star <= ystar ) {
                        p[, ihi] <- p2star
                        y[ihi] <- y2star

                    } else {
                        p[, ihi] <- pstar
                        y[ihi] <- ystar
                    }
                }
            }

            #  Check if YLO improved.
            if ( y[ihi] < ylo ) {
                ylo <- y[ihi]
                ilo <- ihi
            }
            jcount <- jcount - 1
            
            if ( 0 < jcount )
                next

            #  Check to see if minimum reached.
            if ( icount <= kcount ) {
              jcount <- konvge

              x <- sum(y) / dnn
              z <- sum((y - x)^2)
              if ( z <= rq )
                break
            }
        }                               # end inner while loop

        # Factorial tests to check that YNEWLO is a local minimum
        xmin <- p[, ilo]
        ynewlo <- y[ilo]

        if ( kcount < icount ) {
            ifault <- 2
            break
        }

        ifault <- 0

        # Check in all directions with step length
        for (i in 1:n) {
            del <- step[i] * eps
            xmin[i] <- xmin[i] + del
            z <- fn ( xmin )
            icount <- icount + 1
            if ( z < ynewlo ) {
                ifault <- 2
                break
            }
            xmin[i] <- xmin[i] - del - del
            z <- fn ( xmin )
            icount <- icount + 1
            if ( z < ynewlo ) {
                ifault <- 2
                break
            }
            xmin[i] <- xmin[i] + del
        }

        if ( ifault == 0 )
            break

        #  Restart the procedure.
        start <- xmin

        del <- eps
        numres <- numres + 1
    }                                   # end outer while loop

    return(list(xmin = xmin, fmin = ynewlo,
                nfeval = icount, restarts = numres))
} # end of function


.nelminb <- function (x0, fn, lower, upper, tol, maxfeval, step) {
    n <- length(x0)
    if(!is.numeric(lower) || !is.numeric(upper))
        stop("Lower and upper limits must be numeric.", call. = FALSE)
    if (length(lower) == 1) lower <- rep(lower, n)
    if (length(upper) == 1) upper <- rep(upper, n)
    if (length(lower) != n || length(upper) != n)
        stop("Length of 'x0', 'lower', and 'upper' must all be the same.")
    if (any(x0 <= lower) || any (x0 >= upper))
        stop("Starting point must lie in the interior of the bounded region.")
    Trf <- .transfinite(lower, upper, length(x0))
    h <- Trf$h; hinv <- Trf$hinv

    f <- function(x) fn(hinv(x))  # f must be defined on all of R^n
    S <- .nelmin(h(x0), f, tol = tol, maxfeval = maxfeval, step = step)
    S$xmin <- hinv(S$xmin)

    return(S)
}


.transfinite <- function(lower, upper, n = length(lower)) {
    stopifnot(is.numeric(lower), is.numeric(upper))
    if (any(is.na(lower)) || any(is.na(upper)))
        stop("Any 'NA's not allowed in 'lower' or 'upper' bounds.")
    if (length(lower) != length(upper))
        stop("Length of 'lower' and 'upper' bounds must be equal.")
    if (any(lower == upper))
        stop("No component of 'lower' can be equal to the one in 'upper'.")
    if (length(lower) == 1 && n > 1) {
        lower <- rep(lower, n)
        upper <- rep(upper, n)
    } else if (length(lower) != n)
        stop("If 'length(lower)' not equal 'n', then it must be one.")

    low.finite <- is.finite(lower)
    upp.finite <- is.finite(upper)
    c1 <- low.finite & upp.finite    # both lower and upper bounds are finite 
    c2 <- !(low.finite | upp.finite) # both lower and upper bounds infinite
    c3 <- !(c1 | c2) & low.finite    # finite lower bound, infinite upper bound
    c4 <- !(c1 | c2) & upp.finite    # finite upper bound, infinite lower bound

    h <- function(x) {               # h: B --> R^n
        if (any(x < lower) || any (x > upper)) 
            return(rep(NA, n))
        
        hx <- x
        hx[c1] <- atanh(2 * (x[c1] - lower[c1]) / (upper[c1] - lower[c1]) - 1)
        hx[c3] <- log(x[c3] - lower[c3])
        hx[c4] <- log(upper[c4] - x[c4])
        return(hx)
    }

    hinv <- function(x) {            # hinv: R^n --> B
        hix <- x
        hix[c1] <- lower[c1] + (upper[c1] - lower[c1])/2 * (1 + tanh(x[c1]))
        hix[c3] <- lower[c3] + exp(x[c3])
        hix[c4] <- upper[c4] - exp(x[c4])
        return(hix)
    }

    return(list(h = h, hinv = hinv))
}


# nelder_mead <- function(x0, fn, maxfeval = 5000, scale = 1, tol = 1e-10, ...)
# {
#     show <- FALSE
#     if (!is.numeric(x0) || length(x0) < 2)
#         stop("Argument 'x0' must be a numeric vector of length greater 1.")
#     
# 
#     fun <- match.fun(fn)
#     F   <- function(x) scale * fun(x, ...)
# 
#     n <- length(x0)
#     # simplex vertices around x0
#     V <- t(1/(2*n) * cbind(diag(n), rep(-1, n)) + x0)
# 
#     if (show) {
#         P <- Q <- c()
#     }
# 
#     # Function values at vertices
#     Y <- numeric(n+1)
#     for (j in 1:(n+1)) Y[j] <- F(V[j, ])
#     ho <- lo <- which.min(Y)
#     li <- hi <- which.max(Y)
# 
#     for (j in 1:(n+1)) {
#        if (j != lo && j != hi && Y[j] <= Y[li]) li <- j
#        if (j != hi && j != lo && Y[j] >= Y[ho]) ho <- j
#     }
# 
#     cnt <- 0
#     while ( Y[hi] > Y[lo] + tol && cnt < maxfeval ) {
#         S <- numeric(n)
#         for (j in 1:(n+1)) S <- S + V[j,1:n]
#         M <- ( S - V[hi,1:n])/n
#         R <- 2*M - V[hi,1:n]
#         yR <- F(R)
# 
#         if (yR < Y[ho]) {
#            if (Y[li] < yR) {
#               V[hi,1:n] <- R
#               Y[hi] <- yR
#            } else {
#               E <- 2*R - M
#               yE <- F(E)
#               if (yE < Y[li]) {
#                  V[hi,1:n] <- E
#                  Y[hi] <- yE
#               } else {
#                  V[hi,1:n] <- R
#                  Y[hi] <- yR
#               }
#            }
#         } else {
#            if (yR < Y[hi]) {
#               V[hi,1:n] <- R
#               Y[hi] <- yR
#            }
#            C <- (V[hi,1:n] + M)/2
#            yC <- F(C)
#            C2 <- (M + R)/2
#            yC2 <- F(C2)
#            if (yC2 < yC) {
#               C <- C2
#               yC <- yC2
#            }
#            if (yC < Y[hi]) {
#               V[hi,1:n] <- C
#               Y[hi] <- yC
#            } else {
#               for (j in 1:(n+1)) {
#                  if (j != lo) {
#                     V[j,1:n] <- (V[j,1:n] + V[lo,1:n])/2
#                     Z <- V[j,1:n]
#                     Y[j] <- F(Z)
#                  }
#               }
#            }
#         }
# 
#         ho <- lo <- which.min(Y)
#         li <- hi <- which.max(Y)
#         for (j in 1:(n+1)) {
#            if (j != lo && j != hi && Y[j] <= Y[li]) li <- j
#            if (j != hi && j != lo && Y[j] >= Y[ho]) ho <- j
#         }
# 
#         cnt <- cnt + 1
#         if (show) {
#             P <- rbind(P, V[lo, ])
#             Q <- c(Q, Y[lo])
#         }
#     }
# 
#     snorm <- 0
#     for (j in 1:(n+1)) {
#        s <- abs(V[j] - V[lo])
#        if (s >= snorm) snorm <- s
#     }
# 
#     V0 <- V[lo, 1:n]
#     y0 <- Y[lo]
#     dV <- snorm
#     dy <- abs(Y[hi] - Y[lo])
# 
#     if (show) {
#         return(list(xmin = V0, fmin = y0/scale, nfeval = cnt,
#                     dV = dV, dy = dy, P = P, Q = Q))
#     } else {
#         return(list(xmin = V0, fmin = y0/scale, nfeval = cnt))
#     }
# }
