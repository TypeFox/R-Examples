setMethod("getCutoff", "OutlierPCDist", function(obj)
{
    p <- length(getCenter(obj@covobj[[1]]))
    const <- sqrt(qchisq(0.975, p))
   return(rep(const, length(levels(obj@grp))))
})

setMethod("getDistance", "OutlierPCDist", function(obj)
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
OutlierPCDist <- function (x, ...) UseMethod("OutlierPCDist")

OutlierPCDist.formula <- function(formula, data, ..., subset, na.action)
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
    res <- OutlierPCDist.default(x, grouping, ...)

##    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1]] <- as.name("OutlierPCDist")
    res@call <- cl

    res
}


OutlierPCDist.default <- function(x,
                 grouping,
                 control,
                 k,
                 explvar,       ## explained variance for selection of number of components
                 trace=FALSE,
                 ...)
{
    if(is.null(dim(x)))
        stop("x is not a matrix")

    if(missing(control) || is.null(control))
        control = CovControlSest(method="sfast")

    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    grpx <- .getGrouping(grouping, n)

    pca <- PcaClassic(x)
    vars <- getEigenvalues(pca)
    if(missing(k))
        k <- ifelse(missing(explvar), .dimsel(vars, eigen=TRUE), .dimsel.var(vars, eigen=TRUE, explvar=explvar))

    if(k == 1)
        k <- min(4, p)

    data <- as.matrix(getScores(pca)[, 1:k])

    acov <- list()
    wt <- rep(0,n)
    flag <- rep(0,n)
    outliers <- c()
    for(i in 1:grpx$ng)
    {
        class.labels <- which(grpx$grouping == grpx$lev[i])
        class <- data[class.labels,]

        xcov <- CovRobust(class, control=control)

        class.outliers <- which(getFlag(xcov) == 0)
        outliers <- c(outliers, class.labels[class.outliers])
        acov[[i]] <- xcov
        flag[class.labels] <- wt[class.labels] <- getFlag(xcov)
    }

    ret <- new("OutlierPCDist",
                 call=xcall,
                 counts=grpx$counts,
                 grp=grpx$grouping,
                 k=k,
                 covobj=acov,
                 wt = wt,
                 flag=flag,
                 method="PCDIST")

    return (ret)
}

.dimsel.test <- function()
{
    x <- c(10,9,2,1)
    barplot(x)

    ## obvious gap between 2 and 3 eigenvalue
    ##  should select k=2
    .dimsel(x, eigen=TRUE)
}

## Select number of principal components from the screeplot
## If eigen= TRUE, data is a vector containing eigenvalues
.dimsel <- function(data, eigen=FALSE, robust=FALSE)
{
    var.new <- function(x)
    {
        return(if(length(x) == 1) 0 else var(x))
    }

    if(eigen == FALSE)
    {
        if(is.data.frame(data))
            data <- data.matrix(data)
        else if (!is.matrix(data))
            data <- matrix(data, length(data), 1,
              dimnames = list(names(data), deparse(substitute(data))))

        ## drop all rows with missing values (!!) :
        na.x <- !is.finite(data %*% rep(1, ncol(data)))
        ok <- !na.x
        data <- data[ok, , drop = FALSE]
        dx <- dim(data)
        if(!length(dx))
            stop("All observations have missing values!")

        pca <- if(robust) PcaHubert(data) else PcaClassic(data)
        x <- getEigenvalues(pca)
    } else
        x <- data

    n <- length(x)
    profile <- rep(0, n)

    for(i in 1:(n-1))
    {
        x.1 <- x[1:i]
        x.2 <- x[(i+1):n]
        mean.1 <- mean(x.1)
        mean.2 <- mean(x.2)
        var.1 <- var.new(x.1)
        var.2 <- var.new(x.2)
        sd <- ifelse(n <= 2, 1, sqrt(((i-1)*var.1 + (n-i-1)*var.2)/(n-2)))

        profile[i] <- sum(dnorm(x.1,mean.1, sd, log=TRUE)) + sum(dnorm(x.2, mean.2, sd, log=TRUE))
    }
    mean <- mean(x)
    sd <- sd(x)
    profile[n] <- sum(dnorm(x, mean, sd, log=TRUE))
##    screeplot(pca, npcs=100)
    p <- which(profile == max(profile))

##    cat("\n[dimsel] Selected dimension: ", p, "\n")
    p
}

## Select number of principal components which explain at least explvar of variance
## If eigen= TRUE, data is a vector containing eigenvalues
.dimsel.var <- function(data, eigen=FALSE, explvar=0.99)
{
    if(eigen == FALSE)
    {
        if(is.data.frame(data))
            data <- data.matrix(data)
        else if (!is.matrix(data))
            data <- matrix(data, length(data), 1,
              dimnames = list(names(data), deparse(substitute(data))))

        ## drop all rows with missing values (!!) :
        na.x <- !is.finite(data %*% rep(1, ncol(data)))
        ok <- !na.x
        data <- data[ok, , drop = FALSE]
        dx <- dim(data)
        if(!length(dx))
            stop("All observations have missing values!")

        pca <- PcaClassic(data)
        x <- getEigenvalues(pca)
    } else
        x <- data

    x <- x/sum(x)
    csum <- cumsum(x)
    aind <- which(csum >= explvar)
    dimx <- ifelse(length(aind) >= 1, aind[1], csum[length(csum)])
    dimx
}
