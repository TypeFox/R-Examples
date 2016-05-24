

`addrv` <- function (space, FUN = NULL, invars = NULL, name = NULL, ...){
    if (any(class(space) == "ps")) 
        stop("not implemented for class 'ps'")
    if (!is.data.frame(space) | is.null(space$probs)) {
        message("'space' is not a proper probability space")
        stop("see ?probspace")
    }
    bnames <- names(space)[which(names(space) != "probs")]
    out <- subset(space, select = bnames)
    probs <- subset(space, select = probs)
    if (is.null(invars)) 
        invars <- bnames
    if (!is.character(invars)) 
        stop("vars should be a character vector")
    if (!is.null(FUN)) {
        if (is.null(name)) 
            name <- "X"
        temp <- apply(subset(space, select = invars), 1, FUN)
        val <- cbind(out, temp, probs)
        names(val) <- c(bnames, name, "probs")
    }
    else {
        val <- transform(out, ...)
        val$probs <- probs
    }
    return(val)
}



`marginal` <- function (space, vars = NULL){
    if (!is.data.frame(space) | is.null(space$probs)) {
        message("'space' is not a proper probability space")
        stop("see ?probspace")
    }
    if (is.null(vars)) 
        vars <- names(space)[names(space) != "probs"]
    if (!is.character(vars)) {
        stop("'vars' must be a character vector")
    }
    if (length(vars) > 1) {
        res <- aggregate(space$probs, by = as.list(space[, vars]), 
            FUN = sum)
    }
    else {
        res <- aggregate(space$probs, by = list(space[, vars]), 
            FUN = sum)
    }
    names(res) <- c(vars, "probs")
    return(res)
}


`noorder` <- function (space){
    if (!is.data.frame(space)) {
        message("'space' is missing a probs column")
        stop("see ?probspace")
    }
    if (is.null(space$probs)) {
        if (dim(space)[2] < 2) 
            stop("'space' has only one column of outcomes; already unordered")
        n <- names(space)
        res <- unique(data.frame(t(apply(space, 1, sort))))
        names(res) <- n
    }
    else {
        if (dim(space)[2] < 3) 
            stop("'space' has only one column of outcomes; already unordered")
        A <- subset(space, select = -probs)
        probs <- subset(space, select = probs)
        n <- names(A)
        sA <- data.frame(t(apply(A, 1, sort)))
        res <- cbind(sA, probs)
        res <- aggregate(res$probs, by = as.list(sA), sum)
        names(res) <- c(n, "probs")
    }
    return(res)
}
