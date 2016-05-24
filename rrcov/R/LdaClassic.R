##  The S3 version
LdaClassic <- function (x, ...) UseMethod("LdaClassic")

LdaClassic.formula <- function(formula, data, ..., subset, na.action)
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
    res <- LdaClassic.default(x, grouping, ...)

##    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1]] <- as.name("LdaClassic")
    res@call <- cl

##    res$contrasts <- attr(x, "contrasts")
##    res$xlevels <- .getXlevels(Terms, m)
##    res$na.action <- attr(m, "na.action")

    res
}


LdaClassic.default <- function(x,
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

    xcov <- .wcovClass(x, grouping)
    inv <- solve(xcov$wcov)
    ldf <- xcov$means %*% inv
    ldfconst <- diag(log(prior) - ldf %*% t(xcov$means)/2)
    return (new("LdaClassic", call=xcall, prior=prior, counts=counts,
                 center=xcov$means,
                 cov=xcov$wcov,
                 ldf = ldf,
                 ldfconst = ldfconst,
                 method="Linear Discriminant Analysis (LDA)",
                 X=x,
                 grp=g))
}

.wcovClass <- function(x, grouping, method = c("A", "B")){

    xcall <- match.call()
    method <- match.arg(method)
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
    for(i in 1:ng){
        tmpc <- cov.wt(as.matrix(x[which(g == lev[i]),]))
        mX[i,] <- tmpc$center
        covX[,,i] <- tmpc$cov
    }

    if(method == "A"){
        #Method A: pool the covariance matrices
        wcov <- matrix(0,p,p)
        for(i in 1:ng){
            wcov <- wcov + (counts[i]-1)*covX[,,i]
        }
        wcov <- wcov/(n-ng)
    }else if(method == "B"){
        #Method B: center the data and compute the covariance matrix
        #of the centered data
        tmpc <- cov.wt(x - mX[g,])
        tmpc$cov <- tmpc$cov*(n-1)/(n-ng)
        # Add tmpc$center to mX ->mB
        wcov <- tmpc$cov
    }else if(method == "C"){
        stop("Method C not defined for classical estimates")
    }else{
        stop("Unknown method specified: ", method)
    }

    dimnames(wcov) <- list(dimn[[2]], dimn[[2]])
    dimnames(mX) <- list(levels(g), dimn[[2]])

    ans <- list(call=xcall, means=mX, wcov=wcov, method="MLE")
    class(ans) <- "wcov"
    return(ans)
}
