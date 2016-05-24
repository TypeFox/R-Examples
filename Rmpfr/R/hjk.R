####
####  h o o k e j e e v e s . R  Hooke-Jeeves Minimization Algorithm
####

## From: John C Nash <nashjc@uottawa.ca>
## To: Martin Maechler <maechler@stat.math.ethz.ch>, Hans Werner Borchers
## 	<hwborchers@googlemail.com>
## Subject: Re: Rmpfr for optimization? Minor success.
## Date: Tue, 5 Jun 2012 12:37:19 -0400

## Changing to hjk routine was a bit easier to deal with. I found main changes
## were to wrap output with as.numeric() to allow cat() to function.

## I'll get no prizes for tidy code, but it is running an n=2 Chebyquad
## minimization, and seems to be working on an n=6. This may be a good way to
## encourage the sale of cpu power.

## Best, JN

hjkMpfr <- function(par, fn, control = list(), ...) {
    ## Following fails when par is mpfr number JN120605
    ## if (!is.numeric(par))
    ##    stop("Argument 'par' must be a numeric vector.", call. = FALSE)
    n <- length(par)
    if (n == 1)
        stop("For univariate functions use some different method.", call. = FALSE)

    ##-- Control list handling ----------
    cntrl <- list(tol      = 1.e-06,
                  maxfeval = Inf,   # set to Inf if no limit wanted
                  maximize = FALSE, # set to TRUE  for maximization
                  target   = Inf,   # set to Inf for no restriction
                  info     = FALSE) # for printing interim information
    nmsCo <- match.arg(names(control), choices = names(cntrl), several.ok = TRUE)
    if (!is.null(names(control))) cntrl[nmsCo] <- control

    tol      <- cntrl$tol;
    maxfeval <- cntrl$maxfeval
    maximize <- cntrl$maximize
    target   <- cntrl$target
    info <- cntrl$info

    scale <- if (maximize) -1 else 1
    fun <- match.fun(fn)
    f <- function(x) scale * fun(x, ...)

    ##-- Setting steps and stepsize -----
    nsteps <- floor(log2(1/tol))	# number of steps
    steps  <- 2^c(-(0:(nsteps-1)))      # decreasing step size
    dir <- diag(1, n, n)                # orthogonal directions

    x <- par                            # start point
    fx <- f(x)                          # smallest value so far
    fcount <- 1                         # counts number of function calls

    if (info) cat(sprintf("step  nofc %-12s | %20s\n",
                          "fmin", "xpar"))

    ##-- Start the main loop ------------
    ns <- 0
    while (ns < nsteps && fcount < maxfeval && abs(fx) < target) {
        ns <- ns + 1
        hjs    <- .hjsearch(x, f, steps[ns], dir, fcount, maxfeval, target)
        x      <- hjs$x
        fx     <- hjs$fx
        ## found  <- hjs$found
        fcount <- fcount + hjs$finc

	if (info)
	    cat(sprintf("%4d %5d %-12.7g | %-20.15g %-20.15g%s\n",
			ns, fcount, as.numeric(fx/scale),
			as.numeric(x[1]), as.numeric(x[2]),
			if(n > 2)" ...."))
    }

    conv <-
	if (fcount > maxfeval) {
	    warning("Function evaluation limit exceeded -- may not converge.")
	    FALSE
	} else if (abs(fx) > target) {
	    warning("Function exceeds min/max value -- may not converge.")
	    FALSE
	} else
	    TRUE
    fx <- fx / scale                    # undo scaling
    list(par = x, value = fx,
         convergence = conv, feval = fcount, niter = ns)
}

##  Search with a single scale -----------------------------
.hjsearch <- function(xb, f, h, dir, fcount, maxfeval, target) {
    xc <- x <- xb
    finc <- 0
    hje  <- .hjexplore(xb, xc, f, h, dir)
    x    <- hje$x
    fx   <- hje$fx
    found <- hje$found
    finc <- finc + hje$numf

    ## Pattern move
    while (found) {
        d  <- x-xb
        xb <- x
        xc <- x+d
        fb <- fx
        hje  <- .hjexplore(xb, xc, f, h, dir, fb)
        x    <- hje$x
        fx   <- hje$fx
        found <- hje$found
        finc <- finc + hje$numf

        if (!found) {  # pattern move failed
           hje  <- .hjexplore(xb, xb, f, h, dir, fb)
           x    <- hje$x
           fx   <- hje$fx
           found <- hje$found
           finc <- finc + hje$numf
        }
        if (fcount + finc > maxfeval || abs(fx) > target) break
    }

    list(x = x, fx = fx, found=found, finc=finc)
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
    found <- FALSE # do we find a better point ?
    dirh <- h * dir
    fbold <- fx
    for (k in sample.int(n, n)) {       # resample orthogonal directions
        p <- xt + (d. <- dirh[, k])
        fp <- f(p)
        numf <- numf + 1

        if (fp >= fb) {
            p <- xt - d.
            fp <- f(p)
            numf <- numf + 1
        }
        if (fp < fb) {
            found <- TRUE
            xt <- p
            fb <- fp
        }
    }
    if(found) {
        x <- xt
        fx <- fb
    }
    list(x = x, fx = fx, found=found, numf = numf)
}
