setMethod("getCutoff", "OutlierMahdist", function(obj)
{
    p <- length(getCenter(obj@covobj[[1]]))
    const <- sqrt(qchisq(0.975, p))
   return(rep(const, length(levels(obj@grp))))
})

setMethod("getDistance", "OutlierMahdist", function(obj)
{
    dist <- rep(0, length(obj@grp))
    for(i in 1:length(levels(obj@grp)))
    {
        class.labels <- which(obj@grp == levels(obj@grp)[i])
        xcov <- obj@covobj[[i]]
        dist[class.labels] <- sqrt(getDistance(xcov))
    }
    return(dist)
})

##  The S3 version
OutlierMahdist <- function (x, ...) UseMethod("OutlierMahdist")

OutlierMahdist.formula <- function(formula, data, ..., subset, na.action)
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
    res <- OutlierMahdist.default(x, grouping, ...)

##    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1]] <- as.name("OutlierMahdist")
    res@call <- cl

    res
}


OutlierMahdist.default <- function(x,
                 grouping,
                 control,
                 trace=FALSE,
                 ...)
{
    if(is.null(dim(x)))
        stop("x is not a matrix")

    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    if(missing(grouping) || is.null(grouping))
        grouping <- rep(0, n)

    if(length(grouping) == 1) {
        # this is the number of groups and the groups are of equal size
        ng = grouping
        ni = n/ng
        if(ng*ni < n)
            stop("nrow(x) is not divisible by the number of groups")
        grouping <- rep(0,0)
        for(i in 1:ng)
            grouping <- c(grouping, rep(i,ni))
    }else if(length(grouping) > 1 && length(grouping) < n) {
        # grouping contains a vector with the group sizes
        ng <- length(grouping)
        if(sum(grouping) != n)
            stop("nrow(x) is not equal to n1+n2+...+nn")

        gx <- rep(0,0)
        for(i in 1:ng)
            gx <- c(gx, rep(i,grouping[i]))
        grouping <- gx
    }

    if(n != length(grouping))
        stop("nrow(x) and length(grouping) are different")

    g <- as.factor(grouping)
    lev <- lev1 <- levels(g)
    counts <- as.vector(table(g))

    if(any(counts == 0)) {
        warning(paste("group(s)", paste(lev[counts == 0], collapse=" "),"are empty"))
        lev1 <- lev[counts > 0]
        g <- factor(g, levels=lev1)
        counts <- as.vector(table(g))
    }
    proportions <- counts/n
    ng <- length(proportions)
    names(g) <- NULL

    acov <- list()
    wt <- rep(0,n)
    flag <- rep(0,n)
    outliers <- c()
    for(i in 1:length(lev1))
    {
        class.labels <- which(g == lev1[i])
        class <- x[class.labels,]

##        if(nrow(class) < 2*p)
##        {
##            xcov <- CovClassic(class)
##        } else
##        {
##            xcov <- CovRobust(class, control=control)
##        }

        xcov <- CovRobust(class, control=control)

        class.outliers <- which(!getFlag(xcov))
        outliers <- c(outliers, class.labels[class.outliers])
        acov[[i]] <- xcov
        flag[class.labels] <- getFlag(xcov)
        if(!is.null(xcov@wt))
            wt[class.labels] <- xcov@wt
        else
            wt[class.labels] <- getFlag(xcov)
    }

    ret <- new("OutlierMahdist",
                 call=xcall,
                 counts=counts,
                 grp=g,
                 covobj=acov,
                 wt = wt,
                 flag=flag,
                 method=getMeth(xcov))

    return (ret)
}
