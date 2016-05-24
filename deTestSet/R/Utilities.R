### ============================================================================
### first some common functions
### This is the same as in deSolve...
### ============================================================================

## overruling default parameter values 

overrulepar <- function(main, subset, nr) {
 nmsC <- names(main)
    main[(namc <- names(subset))] <- unlist(subset)
    if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
        warning("unknown names in parameter subset: ", paste(noNms, 
            collapse = ", "))
    if (length(main) != nr) 
        stop(paste("length of parameter vector should be", nr))
    return(main)
}

## Checking initial conditions 

checkini <- function (nr, yini, dyini =NULL) {
    if (length(yini) != nr) 
        stop(paste("length of initial state variable vector should be", 
            nr))
    if (!is.null(dyini)) 
        if (length(dyini) != nr) 
            stop(paste("length of initial derivative vector should be", 
                nr))
}

## Overruling integration options
overruleop <- function (main, subset)
{
    nmsC <- names(main)
    main[(namc <- names(subset))] <- subset
    if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
        warning("unknown names in options subset: ", paste(noNms, 
            collapse = ", "))
    return(main)
}

## =============================================================================
## Set the mfrow parameters and whether to "ask" for opening a new device
## =============================================================================

setplotpar <- function(nmdots, dots, nv, ask) {
    if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
      nc <- min(ceiling(sqrt(nv)),3)
      nr <- min(ceiling(nv/nc),3)
      mfrow <- c(nr, nc)
    }
    else if ("mfcol" %in% nmdots)
        mfrow <- rev(dots$mfcol)
    else mfrow <- dots$mfrow

    if (! is.null(mfrow)) {
      mf <- par(mfrow=mfrow)
    }

   ## interactively wait if there are remaining figures
    if (is.null(ask))
      ask <- prod(par("mfrow")) < nv && dev.interactive()

    return(ask)
}

## =============================================================================
## find a variable
## =============================================================================
## uses a loop, such as to keep ordering of inputted vars.

selectvar <- function (which,var) {

    if (!is.numeric(which)) {
        ln <- length(which)
        # keep ordering...
        Select <- NULL
        for ( i in 1:ln) {
          ss <- which(which[i]==var)
          if (length(ss) ==0)
            stop("variable ", which[i], " not in var")
          Select <- c(Select,ss)
        }
    }
    else {
        Select <- which + 1  # "Select now refers to the column number
        if (max(Select) > length(var))
            stop("index in 'which' too large")
        if (min(Select) < 1)
            stop("index in 'which' should be > 0")
    }
  return(Select)
}
