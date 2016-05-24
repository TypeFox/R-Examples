###
### TICTOC.R - Stopwatch timer
###


##-----------------------------------------------------------------------------
tic <- function(gcFirst = FALSE) {
    if (gcFirst == TRUE) {
        gc(verbose = FALSE)
    }
    assign("elapsedTime", proc.time()[3], envir = .pracmaEnv)
    invisible()
}


##-----------------------------------------------------------------------------
toc <- function(echo = TRUE) {
    prevTime <- get("elapsedTime", envir = .pracmaEnv)
    diffTimeSecs <- proc.time()[3] - prevTime
    if (echo) {
        cat(sprintf("elapsed time is %f seconds", diffTimeSecs), "\n")
        return(invisible(diffTimeSecs))
    } else {
        return(diffTimeSecs)
    }
}

