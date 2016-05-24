


`iidspace` <- function (x, ntrials, probs = NULL){
    temp = list()
    for (i in 1:ntrials) {
        temp[[i]] <- x
    }
    res <- expand.grid(temp, KEEP.OUT.ATTRS = FALSE)
    if (is.null(probs)) {
        res$probs <- rep(1/dim(res)[1], dim(res)[1])
    }
    else {
        if (!identical(length(x), length(probs))) {
            stop("'probs' is not the same length as 'outcomes'")
        }
        if (any(probs < 0)) {
            stop("'probs' contains negative values")
        }
        probs <- probs/sum(probs)
        ptemp = list()
        for (i in 1:ntrials) {
            ptemp[[i]] <- probs
        }
        pdf <- expand.grid(ptemp, KEEP.OUT.ATTRS = FALSE)
        pres <- apply(pdf, 1, prod)
        res$probs <- pres
    }
    names(res) <- c(paste(rep("X", ntrials), 1:ntrials, sep = ""), 
        "probs")
    return(res)
}



`is.probspace` <- function (x){
    if (any(class(x) == "ps")) 
        return(TRUE)
    if (!is.data.frame(x) | is.null(x$probs)) 
        return(FALSE)
    if (any(x$probs < 0)) 
        return(FALSE)
    return(TRUE)
}



`probspace` <- function (x, ...)
UseMethod("probspace")



`probspace.default` <- function (x, probs, ...){
    y <- data.frame(x)
    if (missing(probs)) {
        y$probs <- rep(1, dim(y)[1])/dim(y)[1]
    }
    else {
        if (any(probs < 0)) {
            stop("'probs' contains negative values")
        }
        y$probs <- probs/sum(probs)
    }
    return(y)
}


`probspace.list` <- function (x, probs, ...){
    y <- x
    if (missing(probs)) {
        probs <- rep(1, length(y))/length(y)
    }
    else {
        if (any(probs < 0)) {
            stop("'probs' contains negative values")
        }
        probs <- probs/sum(probs)
    }
    res <- list(outcomes = y, probs = probs)
    class(res) <- c("ps", "list")
    return(res)
}
