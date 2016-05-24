##  The S3 version
CSimca <- function (x, ...) UseMethod("CSimca")
CSimca.formula <- function(formula, data=NULL, ..., subset, na.action)
{
    m <- match.call(expand.dots = FALSE)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    if(xint > 0)
        x <- x[, -xint, drop=FALSE]
    res <- CSimca.default(x, grouping, ...)

##    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as 'formula'
    cl <- match.call()
    cl[[1]] <- as.name("CSimca")
    res@call <- cl

##    res$contrasts <- attr(x, "contrasts")
##    res$xlevels <- .getXlevels(Terms, m)
##    res$na.action <- attr(m, "na.action")

    res
}

CSimca.default <- function(x,
                 grouping,
                 prior = proportions,
                 k,
                 kmax=ncol(x),
                 tol = 1.0e-4,
                 trace=FALSE,
                 ...)
{
    if(is.null(dim(x)))
        stop("x is not a matrix")

    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    grpx <- .getGrouping(grouping, n)
    if(!missing(prior))
    {
        if(any(prior < 0) || round(sum(prior), 5) != 1)
            stop("invalid prior")
        if(length(prior) != grpx$ng)
            stop("prior is of incorrect length")
        prior <- prior[grpx$counts > 0]
    }else
        proportions <- grpx$proportions


    if(missing(k))
        k <- rep(0, grpx$ng)
    stopifnot(length(k) == 1 | length(k) == grpx$ng)
    if(length(k) == 1 && k > 0)
        k <- rep(k, grpx$ng)
    else if(length(k) != grpx$ng)
        stop("please provide a value for k")

    acov <- list()
    flag <- rep(0,n)
    for(i in 1:grpx$ng)
    {
        if(trace)
            cat("\nCSimca looping: group=",i,"\n")

        class.labels <- which(grpx$grouping == grpx$lev[i])
        class <- x[class.labels,]

        xpca <- PcaClassic(class, k[i], trace=trace)

        acov[[i]] <- xpca
        flag[class.labels] <- xpca@flag
        k[i] <- xpca@k
    }

    nm <- names(prior)
    if(is.null(nm))
        names(prior) <- grpx$lev

    ret <- new("CSimca",
                 call=xcall,
                 prior=prior,
                 counts=grpx$counts,
                 grp=grpx$grouping,
                 k=k,
                 pcaobj=acov,
                 flag=flag,
                 X=x)

    return (ret)
}
