

`empirical` <- function (x){
    if (any(class(x) == "ps")) 
        stop("not implemented for class 'ps'")
    if (!is.data.frame(x)) 
        message("'x' must be a data frame")
    temp <- x
    n <- dim(temp)[1]
    vars <- names(temp)
    temp$probs <- rep(1, n)/n
    return(marginal(temp))
}



`sim` <- function (x, ...)
UseMethod("sim")


`sim.default` <- function (x, ntrials, ...){
    out <- data.frame(x[, -which(names(x) == "probs")])
    names(out) <- names(x)[-which(names(x) == "probs")]
    p <- x$probs
    d <- dim(x)[1]
    ind <- sample(1:d, size = ntrials, replace = TRUE, prob = p)
    res <- as.data.frame(out[ind, ])
    names(res) <- names(out)
    rownames(res) <- 1:ntrials
    return(res)
}



`sim.ps` <- function (x, ntrials, ...){
    out <- x$outcomes
    p <- x$probs
    d <- length(x$outcomes)
    ind <- sample(1:d, size = ntrials, replace = TRUE, prob = p)
    res <- out[ind]
    return(res)
}
