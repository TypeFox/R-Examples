###
### $Id: tictoc.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Stopwatch timer.
###


##-----------------------------------------------------------------------------
tic <- function(gcFirst = FALSE) {
    if (gcFirst == TRUE) {
        gc(verbose = FALSE)
    }
    assign("savedTime", proc.time()[3], envir = .MatlabNamespaceEnv)
    invisible()
}


##-----------------------------------------------------------------------------
toc <- function(echo = TRUE) {
    prevTime <- get("savedTime", envir = .MatlabNamespaceEnv)
    diffTimeSecs <- proc.time()[3] - prevTime
    if (echo) {
        cat(sprintf("elapsed time is %f seconds", diffTimeSecs), "\n")
        return(invisible())
    } else {
        return(diffTimeSecs)
    }
}

