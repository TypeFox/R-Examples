##
##  h o o k e _ j e e v e s . R  Hooke-Jeeves Minimization
##


hooke_jeeves <- function(x0, f, lb = NULL, ub = NULL, tol = 1e-08,
                          target = Inf, maxfeval = Inf, info = FALSE, ...) {
    if (!is.numeric(x0))
        stop("Argument 'x0' must be a numeric vector.",
              call. = FALSE)
    n <- length(x0)
    if (n == 1)
        stop("Do not use Hooke-Jeeves for univariate functions.",
              call. = FALSE)
    if (!is.function(f) || is.null(f))
        stop("Function argument 'f' must be a valid R function.",
              call. = FALSE)

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    if (is.null(lb) && is.null(ub)) {
        result <- .hj(x0, f, tol, target, maxfeval, info)
    } else {
        result <- .hjb(x0, f, lb, ub, tol, target, maxfeval, info)
    }

    return(result[1:4])
}


.hj <- function(x0, f, tol, target, maxfeval, info) {
    n = length(x0)
    #-- Setting steps and stepsize -----
    nsteps <- floor(log2(1/tol))        # number of steps
    steps  <- 2^c(-(0:(nsteps-1)))      # decreasing step size
    dir <- diag(1, n, n)                # orthogonal directions

    x <- x0                             # start point
    fx <- fbest <- f(x)                 # smallest value so far
    fcount <- 1                         # counts number of function calls

    if (info) cat("step\tnofc\tfmin\txpar\n")   # info header

    #-- Start the main loop ------------
    ns <- 0
    while (ns < nsteps && fcount < maxfeval && abs(fx) < target) {
        ns <- ns + 1
        hjs    <- .hjsearch(x, f, steps[ns], dir, fcount, maxfeval, target)
        x      <- hjs$x
        fx     <- hjs$fx
        sf     <- hjs$sf
        fcount <- fcount + hjs$finc

        if (info)
            cat(ns, "\t",  fcount, "\t", fx, "\t", x[1], "...\n")
    }

    if (fcount > maxfeval) {
        warning("Function evaluation limit exceeded -- may not converge.")
        conv <- 1
    } else if (abs(fx) > target) {
        warning("Function exceeds min/max value -- may not converge.")
        conv <- 1
    } else {
        conv <- 0
    }

    return(list(xmin = x, fmin = fx,
                fcalls = fcount, niter = ns, convergence = conv))
}

##  Search with a single scale -----------------------------
.hjsearch <- function(xb, f, h, dir, fcount, maxfeval, target) {
    x  <- xb
    xc <- x
    sf <- 0
    finc <- 0
    hje  <- .hjexplore(xb, xc, f, h, dir)
    x    <- hje$x
    fx   <- hje$fx
    sf   <- hje$sf
    finc <- finc + hje$numf

    # Pattern move
    while (sf == 1) {
        d  <- x-xb
        xb <- x
        xc <- x+d
        fb <- fx
        hje  <- .hjexplore(xb, xc, f, h, dir, fb)
        x    <- hje$x
        fx   <- hje$fx
        sf   <- hje$sf
        finc <- finc + hje$numf

        if (sf == 0) {  # pattern move failed
           hje  <- .hjexplore(xb, xb, f, h, dir, fb)
           x    <- hje$x
           fx   <- hje$fx
           sf   <- hje$sf
           finc <- finc + hje$numf
        }
        if (fcount + finc > maxfeval || abs(fx) > target) break
    }

    return(list(x = x, fx = fx, sf = sf, finc = finc))
}

##  Exploratory move ---------------------------------------
.hjexplore <- function(xb, xc, f, h, dir, fbold) {
    n <- length(xb)
    x <- xb

    if (missing(fbold)) {
        fb <- f(x)
        numf <- 1
    } else {
        fb <- fbold
        numf <- 0
    }

    fx <- fb
    xt <- xc
    sf <- 0                             # do we find a better point ?
    dirh <- h * dir
    fbold <- fx
    for (k in sample.int(n, n)) {       # resample orthogonal directions
        p1 <- xt + dirh[, k]
        ft1 <- f(p1)
        numf <- numf + 1

        p2 <- xt - dirh[, k]
        ft2 <- f(p2)
        numf <- numf + 1

        if (min(ft1, ft2) < fb) {
            sf <- 1
            if (ft1 < ft2) {
                xt <- p1
                fb <- ft1
            } else {
                xt <- p2
                fb <- ft2
            }
        }
    }

    if (sf == 1) {
        x <- xt
        fx <- fb
    }
    return(list(x = x, fx = fx, sf = sf, numf = numf))
}


.hjb <- function(x0, f, lower, upper, tol, target, maxfeval, info) {
    n <- length(x0)
    # checking lower and upper bounds
    if(!is.numeric(lower) || !is.numeric(upper))
        stop("Lower and upper limits must be numeric.", call. = FALSE)
    if (length(lower) == 1) lower <- rep(lower, n)
    if (length(upper) == 1) upper <- rep(upper, n)
    if (!all(lower <= upper))
        stop("All lower limits must be smaller than upper limits.", call. = FALSE)
    if (!all(lower <= x0) || !all(x0 <= upper))
        stop("Infeasible starting values -- check limits.", call. = FALSE)
    
    #-- Setting steps and stepsize -----
    nsteps <- floor(log2(1/tol))        # number of steps
    steps  <- 2^c(-(0:(nsteps-1)))      # decreasing step size
    dir <- diag(1, n, n)                # orthogonal directions

    x <- x0                            # start point
    fx <- fbest <- f(x)                 # smallest value so far
    fcount <- 1                         # counts number of function calls

    if (info) cat("step\tnofc\tfmin\txpar\n")   # info header

    #-- Start the main loop ------------
    ns <- 0
    while (ns < nsteps && fcount < maxfeval && abs(fx) < target) {
        ns <- ns + 1
        hjs    <- .hjbsearch(x, f, lower, upper,
                            steps[ns], dir, fcount, maxfeval, target)
        x      <- hjs$x
        fx     <- hjs$fx
        sf     <- hjs$sf
        fcount <- fcount + hjs$finc

        if (info)
            cat(ns, "\t",  fcount, "\t", fx, "\t", x[1], "...\n")
    }

    if (fcount > maxfeval) {
        warning("Function evaluation limit exceeded -- may not converge.")
        conv <- 1
    } else if (abs(fx) > target) {
        warning("Function exceeds min/max value -- may not converge.")
        conv <- 1
    } else {
        conv <- 0
    }

    return(list(xmin = x, fmin = fx,
                fcalls = fcount, niter = ns, convergence = conv))
}

##  Search with a single scale -----------------------------
.hjbsearch <- function(xb, f, lo, up, h, dir, fcount, maxfeval, target) {
    x  <- xb
    xc <- x
    sf <- 0
    finc <- 0
    hje  <- .hjbexplore(xb, xc, f, lo, up, h, dir)
    x    <- hje$x
    fx   <- hje$fx
    sf   <- hje$sf
    finc <- finc + hje$numf

    # Pattern move
    while (sf == 1) {
        d  <- x-xb
        xb <- x
        xc <- x+d
        xc <- pmax(pmin(xc, up), lo)
        fb <- fx
        hje  <- .hjbexplore(xb, xc, f, lo, up, h, dir, fb)
        x    <- hje$x
        fx   <- hje$fx
        sf   <- hje$sf
        finc <- finc + hje$numf

        if (sf == 0) {  # pattern move failed
           hje  <- .hjbexplore(xb, xb, f, lo, up, h, dir, fb)
           x    <- hje$x
           fx   <- hje$fx
           sf   <- hje$sf
           finc <- finc + hje$numf
        }
        if (fcount + finc > maxfeval || abs(fx) > target) break
    }

    return(list(x = x, fx = fx, sf = sf, finc = finc))
}

##  Exploratory move ---------------------------------------
.hjbexplore <- function(xb, xc, f, lo, up, h, dir, fbold) {
    n <- length(xb)
    x <- xb

    if (missing(fbold)) {
        fb <- f(x)
        numf <- 1
    } else {
        fb <- fbold
        numf <- 0
    }

    fx <- fb
    xt <- xc
    sf <- 0                             # do we find a better point ?
    dirh <- h * dir
    fbold <- fx
    for (k in sample.int(n, n)) {       # resample orthogonal directions
        p1 <- xt + dirh[, k]
        if ( p1[k] <= up[k] ) {
            ft1 <- f(p1)
            numf <- numf + 1
        } else {
            ft1 <- fb
        }

        p2 <- xt - dirh[, k]
        if ( lo[k] <= p2[k] ) {
            ft2 <- f(p2)
            numf <- numf + 1
        } else {
            ft2 <- fb
        }

        if (min(ft1, ft2) < fb) {
            sf <- 1
            if (ft1 < ft2) {
                xt <- p1
                fb <- ft1
            } else {
                xt <- p2
                fb <- ft2
            }
        }
    }

    if (sf == 1) {
        x  <- xt
        fx <- fb
    }

    return(list(x = x, fx = fx, sf = sf, numf = numf))
}


# hooke_jeeves <- function(x0, f, h = 1, scale = 1,
#                          info = FALSE, tol = .Machine$double.eps^(2/3), ...) {
# 
#     if (!is.numeric(x0))
#         stop("Argument 'x0' must be a numeric vector.", call. = FALSE)
#     x0 <- c(x0)
#     if (!is.numeric(scale) || length(scale) != 1)
#         stop("Argument 'scale' must be a numeric scalar.", call. = FALSE)
#     if (!is.numeric(h) || length(h) != 1 || h <= 0)
#         stop("Argument 'h' must be a positive scalar.", call. = FALSE)
# 
#     fun <- match.fun(f)
#     f <- function(x) scale * fun(x, ...)
# 
#     # Count number of evaluations
#     hjeval <- local({
#         m <- 0
#         function(f, x) {
#             if (missing(f)) return(m)
#             m <<- m + 1
#             return(f(x))
#         }
#     })
# 
#     # Exploration phase
#     hjexplore <- function(x0, f, n, h) {
#         x <- x0; fx <- hjeval(f, x)
#         for (i in sample.int(n, n)) {
#             x[i] <- x[i] + h
#             ff <- hjeval(f, x)
#             if (ff < fx) {
#                 fx <- ff
#             } else {
#                 x[i] <- x0[i] - h
#                 ff <- hjeval(f, x)
#                 if (ff < fx) {
#                     fx <- ff
#                 } else {
#                     x[i] <- x0[i]
#                 }
#             }
#         }
#         return(list(x = x, f = fx))  # shall improvement be returned?
#     }
# 
#     # Main loop
#     n <- length(x0)
#     x <- x0;
#     minf <- hjeval(f, x)
# 
#     # Move phase
#     iter <- 0;
#     while (h > tol) {
#         y <- hjexplore(x, f, n, h)$x
#         if (all(y == x)) {
#             h = h/2
#         } else {
#             x <- 2 * y - x
#             z <- hjexplore(x, f, n, h)$x
#             if (all(z == x)) {
#                 x <- y
#             } else {
#                 x <- z
#             }
#         }
#         iter <- iter + 1
#     }
# 
#     # End program
#     fc <- hjeval(f, x)
#     nc <- hjeval()
#     return(list(xmin = x, fmin = fc, fcalls = nc, niter = iter))
# }
