##  The S3 version
QdaClassic <- function (x, ...) UseMethod("QdaClassic")

QdaClassic.formula <- function(formula, data, ..., subset, na.action)
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
    res <- QdaClassic.default(x, grouping, ...)

##    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as formula
    cl <- match.call()
    cl[[1]] <- as.name("QdaClassic")
    res@call <- cl

##    res$contrasts <- attr(x, "contrasts")
##    res$xlevels <- .getXlevels(Terms, m)
##    res$na.action <- attr(m, "na.action")

    res
}


QdaClassic.default <- function(x,
                 grouping,
                 prior = proportions,
                 tol = 1.0e-4, ...)
{
    if(is.null(dim(x)))
        stop("x is not a matrix")

    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

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

    if(!missing(prior)) {
        if(any(prior < 0) || round(sum(prior), 5) != 1)
            stop("invalid prior")
        if(length(prior) != nlevels(g))
            stop("prior is of incorrect length")
        prior <- prior[counts > 0]
    }
    if(any(counts == 0)) {
        warning(paste("group(s)", paste(lev[counts == 0], collapse=" "),"are empty"))
        lev1 <- lev[counts > 0]
        g <- factor(g, levels=lev1)
        counts <- as.vector(table(g))
    }
    proportions <- counts/n
    ng <- length(proportions)
    names(g) <- NULL
    names(prior) <- levels(g)

    xcov <- .allcovClass(x, grouping)

##    inv <- solve(xcov$wcov)
##    ldf <- xcov$means %*% inv
##    ldfconst <- diag(log(prior) - ldf %*% t(xcov$means)/2)
    return (new("QdaClassic", call=xcall, prior=prior, counts=counts,
                 center=xcov$means,
                 cov=xcov$cov,
                 covinv=xcov$covinv,
                 covdet=xcov$covdet,
                 method="Quadratic Discriminant Analysis (QDA)",
                 control=NULL,
                 X=x,
                 grp=g))
}

.allcovClass <- function(x, grouping){

    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    dimn <- dimnames(x)

    if(!is.factor(g <- grouping))
        g <- as.factor(grouping)
    lev <- levels(g)
    counts <- as.vector(table(g))
    if(any(counts == 0)) {
        stop(paste("group(s)", paste(lev[counts == 0], collapse=" "),"are empty"))
    }
    ng <- length(counts/n)

    # compute group means and covariance matrices for each group
    mX <- matrix(0,ng,p)
    covX <- array(0,c(p,p,ng))
    covInv <- array(0,c(p,p,ng))
    covdet <- vector(mode="numeric", length=ng)
    for(i in 1:ng){
        tmpc <- cov.wt(as.matrix(x[which(g == lev[i]),]))
        mX[i,] <- tmpc$center
        covX[,,i] <- tmpc$cov
        covInv[,,i] <- solve(tmpc$cov)
        covdet[i] <- det(tmpc$cov)
    }

    dimnames(mX) <- list(levels(g), dimn[[2]])
    dimnames(covX) <- list(dimn[[2]], dimn[[2]], levels(g))
    dimnames(covInv) <- list(dimn[[2]], dimn[[2]], levels(g))
    names(covdet) <- levels(g)

    ans <- list(call=xcall, means=mX, cov=covX, covinv=covInv, covdet=covdet)
    class(ans) <- "allcov"
    return(ans)
}
