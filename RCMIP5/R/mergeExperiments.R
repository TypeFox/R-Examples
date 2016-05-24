#' Merge data for two separate experiments
#' 
#' The most common reason to merge experiments is across time periods--e.g.,
#' merging 'historical' with one of the RCPs. This function does that, checking
#' that the merge is appropriate and possible.
#'
#' @param x cmip5data A \code{\link{cmip5data}} object
#' @param y cmip5data A \code{\link{cmip5data}} object
#' @param verbose logical. Print info as we go?
#' @return A \code{\link{cmip5data}} object.
#' @details The variable, units, spatial grid,
#' depths/levels, domain, and model must all match. The timesteps must be
#' identical, and time values non-overlapping. If the time gap between the two experiments
#' is different than their internal timesteps (e.g., if two monthly data objects are
#' separated by more than a month) a warning will be printed.
#' @note This function is 'in the freezer' (not available to users) for now.
#' @note This is an internal RCMIP5 function and not exported.
#' @keywords internal
mergeExperiments <- function(x, y, verbose=FALSE) {
    
    # Sanity checks
    stopifnot(class(x)=="cmip5data" & class(y)=="cmip5data")
    stopifnot(length(verbose)==1 & is.logical(verbose))
    
    if(verbose) cat("Checking that ancillary data are identical\n")
    stopifnot(identical(x$domain, y$domain))
    stopifnot(identical(x$variable, y$variable))
    stopifnot(identical(x$model, y$model))
    stopifnot(identical(x$valUnit, y$valUnit))
    stopifnot(identical(x$lon, y$lon))
    stopifnot(identical(x$lat, y$lat))
    stopifnot(identical(x$Z, y$Z))
    stopifnot(identical(x$lev, y$lev))
    
    # Ensemble check
    if(identical(x$ensembles, y$ensembles)) {
        if(verbose) cat("OK: ensembles match\n")
    } else {
        warning("Ensembles differ between these objects.",
                "\nMerge proceeding but check carefully this is what you want!")   
    }
    
    # Time checks. This is important, and we try to identify obvious problems
    if(verbose) cat("Checking that time data match up\n")
    stopifnot(identical(x$debug$timeFreqStr, y$debug$timeFreqStr))
    
    if(mean(x$time) > mean(y$time)) { # switch them
        temp <- x
        x <- y
        y <- temp
    }
    if(length(intersect(x$time, y$time)) > 0 | max(x$time) > min(y$time)) {
        stop("Overlap between times; can't merge")
    }
    timegap <- min(y$time) - max(x$time)
    tsx <- x$time[2] - x$time[1]
    tsy <- y$time[2] - y$time[1]
    if(isTRUE(all.equal(timegap, tsx)) & isTRUE(all.equal(timegap, tsy))) {
        if(verbose) cat("OK: time gap matches both timesteps\n")
    } else {
        warning("The time gap between objects is ", round(timegap, 3),
                " but their timesteps are ", round(tsx, 3), " and ", round(tsy, 3), 
                "\nMerge proceeding but check carefully this is what you want!")
    }
    
    # Go ahead and merge
    if(verbose) cat("Merging\n")
    x <- addProvenance(x, "Merging with another experiment:")
    x$time <- c(x$time, y$time)
    x$val <- rbind(x$val, y$val)
    x$files <- c(x$files, y$files)
    x$experiment <- paste(x$experiment, y$experiment, sep=".")
    x <- addProvenance(x, y)
    x <- addProvenance(x, "Merge completed")
    x
} # mergeExperiments
