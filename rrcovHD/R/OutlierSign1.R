setMethod("getCutoff", "OutlierSign1", function(obj)
{
    const <- rep(0, length(levels(obj@grp)))
    for(i in 1:length(levels(obj@grp)))
        const[i] <- obj@covobj[[i]]$const
    return(const)
})

setMethod("getDistance", "OutlierSign1", function(obj)
{
    dist <- rep(0, length(obj@grp))
    for(i in 1:length(levels(obj@grp)))
    {
        class.labels <- which(obj@grp == levels(obj@grp)[i])
        dist[class.labels] <- obj@covobj[[i]]$x.dist
    }
    return(dist)
})

setMethod("getCutoff", "OutlierSign2", function(obj)
{
    const <- rep(0, length(levels(obj@grp)))
    for(i in 1:length(levels(obj@grp)))
        const[i] <- obj@covobj[[i]]$const
    return(const)
})

setMethod("getDistance", "OutlierSign2", function(obj)
{
    dist <- rep(0, length(obj@grp))
    for(i in 1:length(levels(obj@grp)))
    {
        class.labels <- which(obj@grp == levels(obj@grp)[i])
        dist[class.labels] <- obj@covobj[[i]]$x.dist
    }
    return(dist)
})

OutlierSign1 <- function (x, ...) UseMethod("OutlierSign1")

OutlierSign1.formula <- function(formula, data, ..., subset, na.action)
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
    res <- OutlierSign1.default(x, grouping, ...)

##    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1]] <- as.name("OutlierSign1")
    res@call <- cl

    res
}


OutlierSign1.default <- function(x,
                 grouping,
                 qcrit = 0.975,
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
    acov <- list()
    wt <- rep(0,n)
    flag <- rep(0,n)
    outliers <- c()
    for(i in 1:grpx$ng)
    {
        class.labels <- which(grpx$grouping == grpx$lev[i])
        class <- x[class.labels,]

        xcov <- .sign1(class, qcrit=qcrit)

        class.outliers <- which(xcov$wfinal01 == 0)
        outliers <- c(outliers, class.labels[class.outliers])
        acov[[i]] <- xcov
        flag[class.labels] <- xcov$wfinal01
        wt[class.labels] <- xcov$wfinal01
    }

    ret <- new("OutlierSign1",
                 call=xcall,
                 counts=grpx$counts,
                 grp=grpx$grouping,
                 covobj=acov,
                 wt = wt,
                 flag=flag,
                 method="SIGN1")

    return (ret)
}

OutlierSign2 <- function (x, ...) UseMethod("OutlierSign2")

OutlierSign2.formula <- function(formula, data, ..., subset, na.action)
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
    res <- OutlierSign2.default(x, grouping, ...)

##    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1]] <- as.name("OutlierSign2")
    res@call <- cl

    res
}


OutlierSign2.default <- function(x,
                 grouping,
                 qcrit = 0.975,
                 explvar = 0.99,
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
    acov <- list()
    wt <- rep(0,n)
    flag <- rep(0,n)
    outliers <- c()
    for(i in 1:grpx$ng)
    {
        class.labels <- which(grpx$grouping == grpx$lev[i])
        class <- x[class.labels,]

        xcov <- .sign2(class, explvar = 0.99, qcrit=qcrit)

        class.outliers <- which(xcov$wfinal01 == 0)
        outliers <- c(outliers, class.labels[class.outliers])
        acov[[i]] <- xcov
        flag[class.labels] <- xcov$wfinal01
        wt[class.labels] <- xcov$wfinal01
    }

    ret <- new("OutlierSign2",
                 call=xcall,
                 counts=grpx$counts,
                 grp=grpx$grouping,
                 covobj=acov,
                 wt = wt,
                 flag=flag,
                 method="SIGN2")

    return (ret)
}

.sign1 <- function (x, qcrit = 0.975)
{
    p = ncol(x)
    n = nrow(x)

    x.mad = apply(x, 2, mad)
    if(any(x.mad == 0))
        stop("More than 50% equal values in one or more variables!")

    x.sc <- scale(x, apply(x, 2, median), x.mad)
    xs <- x.sc/sqrt(apply(x.sc^2, 1, sum))
    xs.evec <- svd(xs)$v
    xs.pc <- x.sc %*% xs.evec
    xs.pcscal <- apply(xs.pc, 2, mad)^2
    xs.pcorder <- order(xs.pcscal, decreasing = TRUE)
    p1 <- min(p-1, n-1)

    dvec <- 1/xs.pcscal[xs.pcorder[1:p1]]
    if(length(dvec) == 1)
        dvec <- as.matrix(dvec)
    covm1 <- xs.evec[, xs.pcorder[1:p1], drop=FALSE] %*%
            diag(dvec) %*%
            t(xs.evec[, xs.pcorder[1:p1], drop=FALSE])
    x.dist <- sqrt(mahalanobis(x.sc, rep(0, p), covm1, inverted = TRUE))
    const <- sqrt(qchisq(qcrit, p1))
    wfinal01 <- (x.dist < const) * 1

    list(wfinal01 = wfinal01,
         x.dist = x.dist,
         const = const)
}

.sign2 <- function (x, explvar=0.99, qcrit=0.975)
{
    p = ncol(x)
    n = nrow(x)

    x.mad = apply(x, 2, mad)
    if (any(x.mad == 0))
        stop("More than 50% equal values in one or more variables!")

    x.sc <- scale(x, apply(x, 2, median), x.mad)
    xs <- x.sc/sqrt(apply(x.sc^2, 1, sum))
    xs.evec <- svd(xs)$v
    xs.pc <- x.sc %*% xs.evec
    xs.pcscal <- apply(xs.pc, 2, mad)^2
    xs.pcorder <- order(xs.pcscal, decreasing = TRUE)

    p1 <- (1:p)[(cumsum(xs.pcscal[xs.pcorder])/sum(xs.pcscal) > explvar)][1]
    x.pc <- x.sc %*% xs.evec[, xs.pcorder[1:p1]]
    xpc.sc <- scale(x.pc, apply(x.pc, 2, median), apply(x.pc, 2, mad))
    xpc.norm <- sqrt(apply(xpc.sc^2, 1, sum))
    xpc.out <- xpc.norm/median(xpc.norm)
    x.dist <- xpc.out * sqrt(qchisq(0.5, p1))
    const <- sqrt(qchisq(qcrit, p1))
    wfinal01 <- rep(0, n)
    wfinal01[x.dist < const] <- 1

    list(wfinal01 = wfinal01,
         x.dist = x.dist,
         const = const)
}
