


`intersect` <- function (x, ...)
UseMethod("intersect")


`intersect.data.frame` <- function (x, y, ...){
    a <- do.call("paste", c(x, sep = "\r"))
    b <- do.call("paste", c(y, sep = "\r"))
    x[match(intersect(a, b), a), ]
}



`intersect.default` <- function (x, y, ...){
    y <- as.vector(y)
    unique(y[match(as.vector(x), y, 0)])
}



`intersect.ps` <- function (x, y, ...){
    a <- do.call("paste", c(x, sep = "\r"))
    b <- do.call("paste", c(y, sep = "\r"))
    e <- match(intersect(a, b), a)
    res <- list(outcomes = x$outcomes[e], probs = x$probs[e])
    class(res) <- c("ps", "list")
    return(res)
}



`setdiff` <- function (x, ...)
UseMethod("setdiff")


`setdiff.data.frame` <- function (x, y, ...){
    a <- do.call("paste", c(x, sep = "\r"))
    b <- do.call("paste", c(y, sep = "\r"))
    x[match(setdiff(a, b), a), ]
}


`setdiff.default` <- function (x, y, ...){
    x <- as.vector(x)
    y <- as.vector(y)
    unique(if (length(x) || length(y)) 
        x[match(x, y, 0) == 0]
    else x)
}


`setdiff.ps` <- function (x, y, ...){
    a <- do.call("paste", c(x, sep = "\r"))
    b <- do.call("paste", c(y, sep = "\r"))
    e <- match(setdiff(a, b), a)
    res <- list(outcomes = x$outcomes[e], probs = x$probs[e])
    class(res) <- c("ps", "list")
    return(res)
}


`subset.ps` <- function (x, subset, ...){
    e <- substitute(subset)
    r <- sapply(x$outcomes, function(t) {
        eval(e, t)
    })
    if (!is.logical(r)) 
        stop("'subset' must be logical")
    res <- list(outcomes = x$outcomes[r & !is.na(r)], probs = x$probs[r & 
        !is.na(r)])
    class(res) <- c("ps", "list")
    return(res)
}


`union` <- function (x, ...)
UseMethod("union")


`union.data.frame` <- function (x, y, ...){
    res <- unique(rbind(as.data.frame(y), x))
    res[order(as.numeric(rownames(res))), ]
}



`union.default` <- function (x, y, ...)
unique(c(as.vector(x), as.vector(y)))


`union.ps` <- function (x, y, ...){
    na <- length(x$outcomes)
    nb <- length(y$outcomes)
    temp <- x
    for (i in 1:nb) {
        temp$outcomes[[na + i]] <- y$outcomes[[i]]
        temp$probs[[na + i]] <- y$probs[[i]]
    }
    a <- do.call("paste", c(temp, sep = "\r"))
    e <- !duplicated(a)
    res <- list(outcomes = temp$outcomes[e], probs = temp$probs[e])
    class(res) <- c("ps", "list")
    return(res)
}
