##
##  h o o k e j e e v e s . R  Hooke-Jeeves Minimization Algorithm
##


hjkb <- function(par, fn, lower = -Inf, upper = Inf, control = list(), ...) {
    if (!is.numeric(par))
        stop("Argument 'par' must be a numeric vector.", call. = FALSE)
    n <- length(par)
    if (n == 1)
        stop("For univariate functions use some different method.", call. = FALSE)

    if(!is.numeric(lower) || !is.numeric(upper))
        stop("Lower and upper limits must be numeric.", call. = FALSE)
    if (length(lower) == 1) lower <- rep(lower, n)
    if (length(upper) == 1) upper <- rep(upper, n)
    if (!all(lower <= upper))
        stop("All lower limits must be smaller than upper limits.", call. = FALSE)
    if (!all(lower <= par) || !all(par <= upper))
        stop("Infeasible starting values -- check limits.", call. = FALSE)
    
   
    #-- Control list handling ----------
    cntrl <- list(tol      = 1.e-06,
                  maxfeval = Inf,       # set to Inf if no limit wanted
                  maximize = FALSE,     # set to TRUE for maximization
                  target   = Inf,       # set to Inf for no restriction
                  info     = FALSE)     # for printing interim information
    nmsCo <- match.arg(names(control), choices = names(cntrl), several.ok = TRUE)
    if (!is.null(names(control))) cntrl[nmsCo] <- control

    tol      <- cntrl$tol;      
    maxfeval <- cntrl$maxfeval
    maximize <- cntrl$maximize
    target   <- cntrl$target
    info     <- cntrl$info

	scale <- if (maximize) -1 else 1
    fun <- match.fun(fn)
    f <- function(x) scale * fun(x, ...)

    #-- Setting steps and stepsize -----
    nsteps <- floor(log2(1/tol))        # number of steps
    steps  <- 2^c(-(0:(nsteps-1)))      # decreasing step size
    dir <- diag(1, n, n)                # orthogonal directions

    x <- par                            # start point
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
            cat(ns, "\t",  fcount, "\t", fx/scale, "\t", x[1], "...\n")
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

    fx <- fx / scale                    # undo scaling
    return(list(par = x, value = fx,
                convergence = conv, feval = fcount, niter = ns))
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
        p <- xt + dirh[, k]
        if (all(p <= up)) {
            ft <- f(p)
            numf <- numf + 1
        } else {
            ft <- fb
        }

        if (ft >= fb) {
            p <- xt - dirh[, k]
            if (all(p >= lo)) {
                ft <- f(p)
                numf <- numf + 1
            }
        }
        if (ft < fb) {
            sf <- 1
            xt <- p
            fb <- ft
        }
    }
    if (sf == 1) {
        x <- xt
        fx <- fb
    }

    return(list(x = x, fx = fx, sf = sf, numf = numf))
}
