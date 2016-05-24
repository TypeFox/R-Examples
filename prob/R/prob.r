
`Prob` <- function (x, ...)
UseMethod("Prob")

`prob` <- function (x, ...){
    message("'prob' is deprecated; use 'Prob' instead.")
    Prob(x, ...)
}

`Prob.default` <- function (x, event = NULL, given = NULL, ...){
    if (is.null(x$probs)) {
        message("'space' is missing a probs column")
        stop("see ?probspace")
    }
    if (missing(event)) {
        r <- TRUE
    }
    else {
        e <- substitute(event)
        r <- eval(e, x, parent.frame())
        if (!is.logical(r)) 
            stop("'event' must evaluate to logical")
        r <- r & !is.na(r)
        if (!isTRUE(all.equal(sum(x$probs), 1))) 
            warning("'space' does not have probability 1.")
    }
    A <- x[r, ]
    if (missing(given)) {
        p <- sum(A$probs)
    }
    else {
        f <- substitute(given)
        g <- eval(f, x, enclos = parent.frame())
        if (!is.logical(g)) {
            if (!is.data.frame(given)) 
                stop("'given' must be data.frame or evaluate to logical")
            B <- given
        }
        else {
            if (missing(event)) 
                stop("'event' must be specified when 'given' is an expression")
            g <- g & !is.na(g)
            B <- x[g, ]
        }
        if (sum(B$probs <= 0)) 
            stop("prob(given) must be positive")
        p <- sum(intersect(A, B)$probs)/sum(B$probs)
    }
    return(p)
}


`Prob.ps` <- function (x, event = NULL, given = NULL, ...){
    if (is.null(x$probs)) {
        message("'space' is missing a probs component")
        stop("see ?probspace")
    }
    if (missing(event)) {
        A <- x
    }
    else {
        e <- substitute(event)
        r <- sapply(x$outcomes, function(t) {
            eval(e, t, enclos=parent.frame())
        })
        if (!is.logical(r)) 
            stop("'event' must evaluate to logical")
        r <- r & !is.na(r)
        if (!isTRUE(all.equal(sum(x$probs), 1))) 
            warning("'space' does not have probability 1.")
        A <- list(outcomes = x$outcomes[r], probs = x$probs[r])
    }
    if (missing(given)) {
        p <- sum(A$probs)
    }
    else {
        f <- substitute(given)
        g <- sapply(x$outcomes, function(t) {
            eval(f, t, enclos=parent.frame())
        })
        if (!is.logical(g)) {
            if (!is.probspace(given)) 
                stop("'given' must be a probspace or evaluate to logical")
            B <- given
        }
        else {
            if (missing(event)) 
                stop("'event' must be specified when 'given' is an expression")
            g <- g & !is.na(g)
            B <- list(outcomes = x$outcomes[g], probs = x$probs[g])
        }
        if (sum(B$probs <= 0)) 
            stop("prob(given) must be positive")
        p <- sum(intersect(A, B)$probs)/sum(B$probs)
    }
    return(p)
}
