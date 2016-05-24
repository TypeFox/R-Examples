##  The S3 version
QdaCov <- function (x, ...) UseMethod("QdaCov")

QdaCov.formula <- function(formula, data, ..., subset, na.action)
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
    res <- QdaCov.default(x, grouping, ...)

##    res$terms <- Terms
    
    ## fix up call to refer to the generic, but leave arg name as 'formula'
    cl <- match.call()
    cl[[1]] <- as.name("QdaCov")
    res@call <- cl

##    res$contrasts <- attr(x, "contrasts")
##    res$xlevels <- .getXlevels(Terms, m)
##    res$na.action <- attr(m, "na.action")

    res
} 

QdaCov.default <- function(x, 
                 grouping, 
                 prior = proportions, 
                 tol = 1.0e-4,
                 method=CovControlMcd(), ...) 
{
    if(is.null(dim(x))) 
        stop("x is not a matrix")

##    method <- match.arg(method)
    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    if(length(grouping) == 1) {
        ## this is the number of groups and the groups are of equal size
        ng = grouping
        ni = n/ng
        if(ng*ni < n)
            stop("nrow(x) is not divisible by the number of groups")
        grouping <- rep(0,0)
        for(i in 1:ng)
            grouping <- c(grouping, rep(i,ni))
    }else if(length(grouping) > 1 && length(grouping) < n) {
        ## grouping contains a vector with the group sizes
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

    if(missing(method))
        method <- "mcd"

    xcov <- .allcovMcd(x, grouping, method, ...)

    return (new("QdaCov", call=xcall, prior=prior, counts=counts, 
                 center=xcov$means, cov=xcov$cov, covinv=xcov$covinv, covdet=xcov$covdet,
                 method=xcov$method, control=xcov$control, X=x, grp=g))
}


.allcovMcd <- function(x, grouping, method=CovControlMcd(), ...){
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
        mcdControl <- .covRobustControl(method)
        mcd <- CovRobust(x[which(g == lev[i]),], control=mcdControl)

        mX[i,] <- getCenter(mcd)
        covX[,,i] <- getCov(mcd)
        covInv[,,i] <- solve(getCov(mcd))
        covdet[i] <- det(getCov(mcd))
    }
    dimnames(mX) <- list(levels(g), dimn[[2]])
    dimnames(covX) <- list(dimn[[2]], dimn[[2]], levels(g))
    dimnames(covInv) <- list(dimn[[2]], dimn[[2]], levels(g))
    names(covdet) <- levels(g)
    ans <- list(call=xcall, means=mX, cov=covX, covinv=covInv, covdet=covdet, method=getMeth(mcd), control=mcdControl)
    class(ans) <- "allcov"
    return(ans)
}
