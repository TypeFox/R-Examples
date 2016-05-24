nm <- function (x, ...) 
    UseMethod("nm")

nm.default <- function(x, grouping, gamma=0, ...)
{
    xlist <- split(data.frame(x),as.factor(grouping))
    x <- t(sapply(xlist, colMeans))
    grouping <- as.factor(levels(grouping))
    k <- 1
    if(gamma != 0) k <- length(grouping)
    if(nrow(x) == 1) x <- matrix(x, ncol=1)
    res <- sknn(x, grouping, k=k, gamma=gamma, ...)
    return(res)
}

### nm bei verschiedenen Eingabeformaten:
nm.formula <- function (formula, data = NULL, ..., subset, na.action = na.fail) 
{
    m <- match.call(expand.dots = FALSE)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0) 
        x <- x[, -xint, drop = FALSE]
    res <- nm.default(x, grouping, ...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1]] <- as.name("nm")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- .getXlevels(Terms, m)
    res$na.action <- attr(m, "na.action")
    res
}

nm.matrix<-function (x, grouping, ..., subset, na.action = na.fail) 
{
    if (!missing(subset)) {
        x <- x[subset, , drop = FALSE]
        grouping <- grouping[subset]
    }
    if (!missing(na.action)) {
        dfr <- na.action(structure(list(g = grouping, x = x), 
            class = "data.frame"))
        grouping <- dfr$g
        x <- dfr$x
    }
    res <- nm.default(x, grouping, ...)
    cl <- match.call()
    cl[[1]] <- as.name("nm")
    res$call <- cl
    res
}



nm.data.frame<-function (x, ...) 
{
   res <- nm.matrix(structure(data.matrix(x), class = "matrix"), 
        ...)
    cl <- match.call()
    cl[[1]] <- as.name("nm")
    res$call <- cl
    res
}
