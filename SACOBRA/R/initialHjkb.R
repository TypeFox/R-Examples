##
## Hooke & Jeeves for initial phase of COBRA
##  h o o k e j e e v e s . R  Hooke-Jeeves Minimization Algorithm
##
initialHjkb <- function(par, fn, lower = -Inf, upper = Inf, control = list(), ...) {
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
    y <- f(x)
	  violatedLines = which(y[2:length(y)]>0) + 1
  	penalty = length(violatedLines)*y[1] + y[violatedLines]*y[1]
	  y[1] = y[1] + sum(penalty) 
    fx <- fbest <- y[1]                 # smallest value so far
    fcount <- 1                         # counts number of function calls
  
    #PK new structures
    xArchive = data.frame()
	  yArchive = c()
	  constraintArchive = data.frame()
    importances = data.frame()
    historyData = list(xArchive=xArchive, yArchive=yArchive, 
                       constraintArchive=constraintArchive, 
                       fEvals=fcount, bestX=x, bestY=y)

    if (info) cat("step\tnofc\tfmin\txpar\n")   # info header

    #-- Start the main loop ------------
    ns <- 0
    while (ns < nsteps && fcount < maxfeval && abs(fx) < target) {
        ns <- ns + 1
        hjs    <- initialHjbsearch(x, f, lower, upper,
                            steps[ns], dir, fcount, maxfeval, target, historyData)
        x      <- hjs$x
        fx     <- hjs$fx
        sf     <- hjs$sf
        fcount <- fcount + hjs$finc
        historyData <- hjs$historyData

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
                convergence = conv, feval = fcount, niter = ns, historyData=historyData))
}

##  Search with a single scale -----------------------------
initialHjbsearch <- function(xb, f, lo, up, h, dir, fcount, maxfeval, target, historyData) {
    x  <- xb
    xc <- x
    sf <- 0
    finc <- 0
    hje  <- initialHjbexplore(xb, xc, f, lo, up, h, dir, historyData=historyData)
    x    <- hje$x
    fx   <- hje$fx
    sf   <- hje$sf
    finc <- finc + hje$numf
    historyData <- hje$historyData

    # Pattern move
    while (sf == 1) {
        d  <- x-xb
        xb <- x
        xc <- x+d
        fb <- fx
        hje  <- initialHjbexplore(xb, xc, f, lo, up, h, dir, fb, historyData)
        x    <- hje$x
        fx   <- hje$fx
        sf   <- hje$sf
        finc <- finc + hje$numf
        historyData <- hje$historyData

        if (sf == 0) {  # pattern move failed
           hje  <- initialHjbexplore(xb, xb, f, lo, up, h, dir, fb, historyData)
           x    <- hje$x
           fx   <- hje$fx
           sf   <- hje$sf
           finc <- finc + hje$numf
           historyData <- hje$historyData
        }
        if (fcount + finc > maxfeval || abs(fx) > target) break
    }

    return(list(x = x, fx = fx, sf = sf, finc = finc, historyData=historyData))
}

##  Exploratory move ---------------------------------------
initialHjbexplore <- function(xb, xc, f, lo, up, h, dir, fbold, historyData) {
    n <- length(xb)
    x <- xb    
    y = Inf

    if (missing(fbold)) {
        y <- f(x)
        historyData$yArchive <- c(historyData$yArchive, y[1])
        violatedLines = which(y[2:length(y)]>0) + 1
        penalty = length(violatedLines)*y[1] + y[violatedLines]*y[1]
        y[1] = y[1] + sum(penalty)
        historyData$xArchive <- rbind(historyData$xArchive, x)
        historyData$constraintArchive <- rbind(historyData$constraintArchive, y[2:length(y)])
        fb <- y[1]
        numf <- 1
        historyData$fEvals <- historyData$fEvals + 1
    } else {
        fb <- fbold
        y <- fbold
        numf <- 0
    }

    fx <- fb
    xt <- xc
    sf <- 0                             # do we find a better point ?
    dirh <- h * dir
    fbold <- fx

    # OPTIMIZATION LOOP
    for (k in sample.int(n, n)) {       # resample orthogonal directions
        p <- xt + dirh[, k]
        if (all(p <= up)) {
            y <- c()
            y <- f(p)
            violatedLines = which(y[2:length(y)]>0) + 1
            penalty = length(violatedLines)*y[1] + y[violatedLines]*y[1]
            historyData$yArchive <- c(historyData$yArchive, y[1])
            y[1] = y[1] + sum(penalty)
            names(p) = paste("x",1:length(p),sep="")
            numf <- numf + 1
            historyData$fEvals <- historyData$fEvals + 1
            colnames(historyData$xArchive) <- paste("x",1:length(p),sep="")
            historyData$xArchive <- rbind(historyData$xArchive, t(p))
            historyData$constraintArchive <- rbind(historyData$constraintArchive, y[2:length(y)])
            ft <- y[1]
            
        } else {
            ft <- fb
        }

        if (ft >= fb) {
            p <- xt - dirh[, k]
            if (all(p >= lo)) {
                y <- c()
                y <- f(p)
                violatedLines = which(y[2:length(y)]>0) + 1
                penalty = length(violatedLines)*y[1] + y[violatedLines]*y[1]
                historyData$yArchive <- c(historyData$yArchive, y[1])
                y[1] = y[1] + sum(penalty) # penalty function
                names(p) = paste("x",1:length(p),sep="")
                numf <- numf + 1
                historyData$fEvals <- historyData$fEvals + 1
                colnames(historyData$xArchive) <- paste("x",1:length(p),sep="")
                historyData$xArchive <- rbind(historyData$xArchive, t(p))
                historyData$constraintArchive <- rbind(historyData$constraintArchive, y[2:length(y)])            
                ft <- y[1]
            }
        }
        if (ft < fb) {
            sf <- 1
            xt <- p
            fb <- ft
            historyData$bestX=p
            historyData$bestY=fb
        }
    }
    if (sf == 1) {
        x <- xt
        fx <- fb
    }

    return(list(x = x, fx = fx, sf = sf, numf = numf, historyData=historyData))
}
