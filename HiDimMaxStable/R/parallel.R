
check.parallel<-function(parallel=TRUE) {
    if(parallel) {
        if (requireNamespace("snowfall", quietly = TRUE)) {
            if(!snowfall::sfIsRunning()) stop("You must call sfInit first")
        } else {
            stop("Package snowfall must be installed.")
        }
    }
}
