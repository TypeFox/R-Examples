multispecies <- 
function(..., FUN=list(mean=mean, sd=sd), na.rm=FALSE) {
    
    dt <- list(...)
    
    # Some error checking
    nvars <- length(dt)
    test <- 1
    for (i in 1:nvars) test <- test * is.numeric(dt[[i]])
    if (!as.logical(test)) stop("Input data must be numeric.")

    if (nvars == 1) {
        if (!(is.data.frame(dt[[1]]) | is.matrix(dt[[1]]))) {
            stop("Single input must be a numeric data.frame or matrix.")
        } else {
            ndt <- dt[[1]]
        }
    } else {
        ndt <- matrix(NA, length(dt[[1]]), nvars)
        for (i in 1:nvars) {
            if (!is.vector(dt[[i]])) {
                stop(paste("Using a '", class(dt[[i]]), 
                           "' with multiple input arguments. Should use ",
                           "only 'vectors' or a single argument of class '",
                           class(dt[[i]]), "' with all data.", sep=''))
            } else {
                ndt[,i] <- dt[[i]]
            }
        }
    }

    result <- matrix(NA, nrow(ndt), length(FUN))
    for (i in 1:length(FUN))
        result[,i] <- apply(ndt, 1, FUN[[i]], na.rm=na.rm)
    if (!is.null(names(FUN))) colnames(result) <- names(FUN)

    result
}
