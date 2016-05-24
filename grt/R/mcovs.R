mcovs <- function(x, ...) UseMethod("mcovs")

mcovs.formula <- function(formula, data, pooled=TRUE, ...)
{
    m <- match.call(expand.dots = FALSE)
    m$... <- NULL
    m$pooled <- NULL
    m[[1L]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xint <- match("(Intercept)", colnames(x), nomatch=0L)
    if(xint > 0) x <- x[, -xint, drop=FALSE]
    res <- mcovs.default(x = x, grouping = grouping, pooled = pooled, ...)
    #res$terms <- Terms
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1L]] <- as.name("mcovs")
    #res$call <- cl
    #res$contrasts <- attr(x, "contrasts")
    #res$xlevels <- .getXlevels(Terms, m)
    res
}

mcovs.default <- function(x, grouping, pooled=TRUE, ...){
    if(!is.matrix(x)) x <- as.matrix(x)
    if(any(!is.finite(x))) stop("infinite, NA or NaN values in 'x'")
    n <- nrow(x)
    if(n != length(grouping))
        stop("nrow(x) and length(grouping) are different")
    p <- ncol(x)
    g <- as.factor(grouping)
    counts <- as.vector(table(g))
    k <- length(counts)
    x <- as.data.frame(x)
    gmeans <- lapply(split(x, grouping), colMeans)
    covList <- if(pooled == FALSE) {
        lapply(split(x, grouping), cov)
    } else {
        list( pooled = matrix(
            (sapply(split(x, grouping), cov) %*% counts)/(n - k), 
            ncol=p)
        )
    }
    res <- list(N = n, counts = counts, lev = levels(g), 
                means = gmeans, covs=covList)
    res
}

