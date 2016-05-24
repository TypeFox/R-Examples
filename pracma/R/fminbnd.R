##
##  f m i n b n d . R  Brent's Minimization Algorithm
##


fminbnd <- function(f, a, b, maxiter = 1000, maximum = FALSE,
                    tol = 1e-07, rel.tol = tol, abs.tol = 1e-15, ...) {
    stopifnot(is.numeric(a), length(a) == 1,
              is.numeric(b), length(b) == 1)
    if (a >= b)
        stop("Interval end points must fulfill  a < b !")
    
    fun <- match.fun(f)
    if (maximum)
        f <- function(x) -fun(x, ...)
    else
        f <- function(x)  fun(x, ...)
    
    # phi is the square of the inverse of the golden ratio.
    phi <- 0.5 * ( 3.0 - sqrt ( 5.0 ) )

    xl <- a; xu <- b
    xmin <- xl + phi*(xu-xl)
    fmin <- f(xmin)
    nfeval   <- 1
    
    step <- step2 <- 0.0

    xmin1 <- xmin; fmin1 <- fmin
    xmin2 <- xmin; fmin2 <- fmin

    converged <- FALSE
    iter <- 0
    while (iter < maxiter) {
        pp <- qq <- 0.0

        tolx <- rel.tol * abs(xmin) + abs.tol

        xm = (xu+xl)/2

        if (abs(xmin - xm) <= 2*tolx - (xu-xl)/2) {
            converged <- TRUE
            break
        }

        iter <- iter + 1
        
        if (abs(step2) > tolx) {
            rr <- (xmin - xmin1) * (fmin - fmin2)
            qq <- (xmin - xmin2) * (fmin - fmin1)
            pp <- (xmin - xmin2) * qq - (xmin - xmin1) * rr
            qq <- 2*(qq - rr)

            if (qq > 0) {
                pp <- -pp
            } else {
                qq <- -qq
            }
        }

        if (abs(pp) < abs(qq*step2/2) && pp < qq*(xu-xmin) && pp < qq*(xmin-xl)) {
            step2 <- step
            step <- pp/qq

            xtemp <- xmin + step
            if ((xtemp - xl) < 2*tolx || (xu - xtemp) < 2*tolx) {
                step <- if (xmin < xm) tolx else -tolx
            }
        } else {
            step2 <- if (xmin < xm) xu - xmin else xl - xmin
            step <- phi * step2
        }

        if (abs(step) >= tolx) {
            xnew <- xmin + step
        } else {
            xnew <- xmin + (if (step > 0) tolx else -tolx)
        }

        fnew <- f(xnew)
        nfeval <- nfeval + 1

        if (fnew <= fmin) {
            if (xnew < xmin) {
                xu <- xmin
            } else {
                xl <- xmin
            }
            xmin2 <- xmin1; fmin2 <- fmin1
            xmin1 <- xmin;  fmin1 <- fmin
            xmin  <- xnew;  fmin <- fnew
        } else {
            if (xnew < xmin) {
                xl <- xnew
            } else {
                xu <- xnew
            }
            if (xnew <= xmin1 || xmin1 == xmin) {
                xmin2 = xmin1; fmin2 = fmin1
                xmin1 = xnew;  fmin1 = fnew
            } else if (xnew <= xmin2 || xmin2 == xmin || xmin2 == xmin1) {
                xmin2 <- xnew; fmin2 <- fnew
            }
        }
    }
    return(list(xmin = xmin, fmin = fmin,
                niter = iter, estim.prec = abs(step)))
}


# fminbnd <- function(f, a, b, ..., maxiter = 1000, maximum = FALSE,
#                     tol = .Machine$double.eps^(2/3)) {
#     stopifnot(is.numeric(a), length(a) == 1,
#               is.numeric(b), length(b) == 1)
#     if (a >= b)
#         stop("Interval end points must fulfill  a < b !")
# 
#     fun <- match.fun(f)
#     if (maximum)
#         f <- function(x) -fun(x, ...)
#     else
#         f <- function(x)  fun(x, ...)
# 
#     # phi is the square of the inverse of the golden ratio.
#     phi <- 0.5 * ( 3.0 - sqrt ( 5.0 ) )
# 
#     # Set tolerances
#     tol1 <- 1 + eps()
#     eps0 <- sqrt(eps())
#     tol3 <- tol / 3
# 
#     sa <- a; sb <- b
#     x  <- sa + phi * ( b - a )
#     fx <- f(x)
#     v  <- w  <- x
#     fv <- fw <- fx
#     e  <- 0.0;
# 
#     niter <- 1
#     while ( niter <= maxiter ) {
#         xm <-  0.5 * ( sa + sb )
#         t1 <- eps0 * abs ( x ) + tol/3
#         t2 <-  2.0 * t1
# 
#         #  Check the stopping criterion.
#         if ( abs ( x - xm ) <= t2 - (dx <- ( sb - sa ) / 2 ) ) break
# 
#         r <- 0.0
#         p <- q <- r
# 
#         # Fit a parabola.
#         if ( t1 < abs ( e ) ) {
#             r <- ( x - w ) * ( fx - fv )
#             q <- ( x - v ) * ( fx - fw )
#             p <- ( x - v ) * q - ( x - w ) * r
#             q <- 2.0 * ( q - r );
#             
#             if ( 0.0 < q ) p <- - p
#             
#             q <- abs ( q )
#             r <- e
#             e <- d
#         }
# 
#         # Is the parabola acceptable
#         if ( abs ( p ) < abs ( 0.5 * q * r ) &&
#                  q * ( sa - x ) < p          &&
#                  p < q * ( sb - x ) ) {
#              #  Take the parabolic interpolation step.
#              d <- p / q
#              u <- x + d
# 
#              #  F must not be evaluated too close to a or b.
#              if ( ( u - sa ) < t2 | ( sb - u ) < t2 ) {
#                  d <- if (x < xm) t1 else -t1
#              }
# 
#          } else {
#          #  A golden-section step.
#             e <- if (x < xm) sb - x else a - x
#             d <- phi * e
#         }
# 
#         #  F must not be evaluated too close to X.
#         if ( t1 <= abs ( d ) ) {
#             u = x + d
#         } else if ( 0.0 < d ) {
#             u = x + t1
#         } else {
#             u = x - t1
#         }
# 
#         fu = f ( u )
# 
#         #  Update a, b, v, x, and x.
#         if ( fu <= fx ) {
#             if ( u < x ) sb <- x
#             else         sa <- x
#         
#             v <- w; fv <- fw
#             w <- x; fw <- fx
#             x <- u; fx <- fu
# 
#         } else {
#             if ( u < x ) sa <- u
#             else         sb <- u
# 
#             if ( fu <= fw || w == x ) {
#                 v <- w; fv <- fw
#                 w <- u; fw <- fu
#             } else if ( fu <= fv || v == x || v== w ) {
#                 v <- u; fv <- fu
#             }
#         }
#         niter <- niter + 1
# 
#     } #endwhile
# 
#     if (niter > maxiter)
#         warning("No. of max. iterations exceeded; no convergence reached.")
# 
#     if (maximum) fx <- -fx
#     return( list(xmin = x, fmin = fx, niter = niter, estim.prec = dx) )
# }
