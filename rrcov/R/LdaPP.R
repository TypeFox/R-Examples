##  The S3 version
LdaPP <- function (x, ...) UseMethod("LdaPP")

LdaPP.formula <- function(formula, data, subset, na.action, ...)
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
    res <- LdaPP.default(x, grouping, ...)

##    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as 'formula'
    cl <- match.call()
    cl[[1]] <- as.name("LdaPP")
    res@call <- cl

##    res$contrasts <- attr(x, "contrasts")
##    res$xlevels <- .getXlevels(Terms, m)
##    res$na.action <- attr(m, "na.action")

    res
}

LdaPP.default <- function(x,
                 grouping,
                 prior = proportions,
                 tol = 1.0e-4,
                 method = c("huber", "mad", "sest", "class"),
                 optim = FALSE,
                 trace=FALSE, ...)
{
    if(is.null(dim(x)))
        stop("x is not a matrix")

    method <- match.arg(method)
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
    #############################################################

    mrob <- if(method == "huber") .mrob.huber
            else if(method == "mad") .mrob.mad
            else if(method == "sest") .mrob.s
            else if(method == "class") .mrob.sd

    ## We have only two groups!
    n1 <- counts[1]
    n2 <- counts[2]
##    a1 <- (n1-1)/(n-2)
    a1 <- n1/n
    a2 <- 1 - a1
    raw.ldf <- ldf <- matrix(0, nrow=2, ncol=p)
    raw.ldfconst <- ldfconst <- rep(0,2)

    alpha <- .projpp(x, grp=g, a1, a2, mrob)
    dx <- .det.a0(alpha, x, grp=g, a1=a1, a2=a2, mrob=mrob, prior=prior)
    raw.ldf[1,] <- dx$alpha
    raw.ldfconst[1] <- dx$alpha0
    ldf <- raw.ldf
    ldfconst <- raw.ldfconst

    if(optim)
    {
        sol.1 <- optim(dx$alpha, .auxPP, x=x, grp=g, a1=a1, a2=a2, mrob=mrob)
        if(sol.1$convergence != 0)
            sol.1 <- optim(sol.1$par, .auxPP, x=x, grp=g, a1=a1, a2=a2, mrob=mrob)
        dx <- .det.a0(sol.1$par, x, grp=g, a1=a1, a2=a2, mrob=mrob, prior=prior)
        ldf[1,] <- dx$alpha
        ldfconst[1] <- dx$alpha0
    }

    return (new("LdaPP", call=xcall, prior=prior, counts=counts,
                 raw.ldf=raw.ldf, raw.ldfconst=raw.ldfconst,
                 ldf=ldf, ldfconst=ldfconst,
                 method=method, X=x, grp=g))
}

.projpp <- function(x, grp, a1, a2, mrob)
{
    lev <- levels(grp)
    counts <- as.vector(table(grp))
    p <- ncol(x)
    n1 <- counts[1]
    n2 <- counts[2]

    X1  <- x[grp == lev[1],]
    X2  <- x[grp == lev[2],]

    inull <- 0
    dmax <- 0
    alpha <- NULL
    imax <- 0
    jmax <- 0
    for(i in 1:n1)
    {
        for(j in 1:n2)
        {
            aux <- sqrt(sum((X1[i,] - X2[j,])**2))
            dx  <- (X1[i,] - X2[j,])/aux
            px1  <- as.vector(t(X1%*%dx))
            px2  <- as.vector(t(X2%*%dx))
            m1 <- mrob(px1)
            m2 <- mrob(px2)
            if(is.null(m1$mu) | is.null(m1$s) | is.null(m2$mu) | is.null(m2$s))
            {
                inull <- inull + 1
                if(inull > 10)
                    stop("Too many unsuccessful S-estimations!")
            }else
            {
                xdist <- (m1$mu - m2$mu)**2/(m1$s**2*a1 + m2$s**2*a2)
                if(xdist > dmax)
                {
                    dmax <- xdist
                    alpha <- dx
                    imax <- i
                    jmax <- j
                }
            }
        }
    }
    alpha
}

.det.a0 <- function(alpha, x, grp, a1, a2, mrob, prior)
{
    lev <- levels(grp)
    counts <- as.vector(table(grp))
    p <- ncol(x)
    n1 <- counts[1]
    n2 <- counts[2]
    n <- n1 + n2

    alpha  <- alpha/sqrt(sum(alpha**2))
    px <- as.vector(t(x%*%alpha))
    m1 <- mrob(px[grp == lev[1]])
    m2 <- mrob(px[grp == lev[2]])
    if(m1$mu < m2$mu)
    {
        alpha <- -1*alpha
        m1$mu <- -1*m1$mu
        m2$mu <- -1*m2$mu
    }

    alpha0 <- -((m1$mu + m2$mu)/2 + log(prior[2]/prior[1])/(m1$mu - m2$mu)*(a1*m1$s**2 + a2*m2$s**2))
    names(alpha0) <- NULL
    list(alpha=alpha, alpha0=alpha0)
}

##
## A function to be maximized by Nelder and Mead algorithm in optim(),
##  with first argument the vector of parameters over which
##  maximization is to take place. It returns a scalar result.
##
##  The initial values are from the projection algorithm
##
.auxPP <- function(alpha, x, grp, a1, a2, mrob)
{
    lev <- levels(grp)
    counts <- as.vector(table(grp))
    p <- ncol(x)
    n1 <- counts[1]
    n2 <- counts[2]
    n <- n1 + n2

    px <- as.vector(t(x%*%alpha))
    m1 <- mrob(px[grp == lev[1]])
    m2 <- mrob(px[grp == lev[2]])
    ret <- -(m1$mu - m2$mu)**2/(m1$s**2*a1 + m2$s**2*a2)

    if(is.null(ret))
        ret <- NA
    if(length(ret) == 0)
    {
        cat("\nERROR: length of distance for 'optim' equal to 0: ", ret, "\n")
        ret <- NA
    }

    ret
}


.mrob.sd <- function(y)
{
    list(mu=mean(y, na.rm=TRUE), s=sd(y, na.rm=TRUE))
}

.mrob.mad <- function(y)
{
    list(mu=median(y, na.rm=TRUE), s=mad(y, na.rm=TRUE))
}

.mrob.s <- function(y)
{
    aux <- lmrob(y ~ 1)
    mu  <- aux$init.S$coef

    names(mu) <- NULL
    s   <- aux$init.S$scale

##    mm=CovSest(as.matrix(y))
##    list(mu=mm@center, s=mm@cov[1,1])

    list(mu=mu, s=s)
}

.mrob.huber <- function(y, k1=1.5, k2=2, tol=1e-06, iter=8)
{
    if(any(i <- is.na(y)))
        y <- y[!i]
    n <- length(y)

    ## replicate the huberM() function from robustbase and add an
    ##  iter parameter to be able to reproduce Ana's results
#    mx1 <- huberM(y, k=k1, tol=tol)
    mx1 <- .huberM(y, k=k1, tol=tol, iter=iter)

    k <- k2
    mu <- median(y)
    s0 <- mx1$s
    th <- 2 * pnorm(k) - 1
    beta <- th + k^2 * (1 - th) - 2 * k * dnorm(k)
    it <- 0:0
    repeat
    {
        it <- it + 1:1
        yy <- pmin(pmax(mu - k * s0, y), mu + k * s0)
        ss <- sum((yy - mu)^2)/n
        s1 <- sqrt(ss/beta)
##        if(abs(s0-s1) < tol || (iter > 0 & it > iter))
        if(it > iter)
            break
        s0 <- s1
    }
    list(mu=mx1$mu, s=s0, it1=mx1$it, it2=it)
}

##
## This is from robustbase - I added 'iter' parameter,
##  only to be able to reproduce exactly na Pires' LdaPP
##  FIXME: remove later
##
.huberM <- function (x, k = 1.5, weights = NULL, tol = 1e-06, mu = if (is.null(weights)) median(x) else wgt.himedian(x,
    weights), s = if (is.null(weights)) mad(x, center = mu) else wgt.himedian(abs(x -
    mu), weights), warn0scale = getOption("verbose"), iter=8)
{
    if (any(i <- is.na(x))) {
        x <- x[!i]
        if (!is.null(weights))
            weights <- weights[!i]
    }
    n <- length(x)
    sum.w <- if (!is.null(weights)) {
        stopifnot(is.numeric(weights), weights >= 0, length(weights) ==
            n)
        sum(weights)
    }
    else n
    it <- 0:0
    if (sum.w == 0)
        return(list(mu = NA, s = NA, it = it))
    if (s <= 0) {
        if (s < 0)
            stop("negative scale 's'")
        if (warn0scale && n > 1)
            warning("scale 's' is zero -- returning initial 'mu'")
    }
    else {
        wsum <- if (is.null(weights))
            sum
        else function(u) sum(u * weights)
        repeat {
            it <- it + 1:1
            y <- pmin(pmax(mu - k * s, x), mu + k * s)
            mu1 <- wsum(y)/sum.w
##            if (abs(mu - mu1) < tol * s)
            if (it > iter)
                break

            mu <- mu1
        }
    }
    list(mu = mu, s = s, it = it)
}

##
## Predict method for LdaPP - additional parameter raw=FALSE.
##  If set to TRUE, the prediction will be done using the raw estimates
##  (obtained by the first approximation algorithm).
##
setMethod("predict", "LdaPP", function(object, newdata, raw=FALSE){

    ct <- FALSE
    if(missing(newdata))
    {
        newdata <- object@X         # use the training sample
        ct <- TRUE                  # perform cross-validation
    }

    x <- as.matrix(newdata)

    if(length(object@center)>0 & ncol(x) != ncol(object@center) | ncol(x) != ncol(object@ldf))
        stop("wrong number of variables")

    ldf <- if(raw) object@raw.ldf else object@ldf
    ldfconst <- if(raw) object@raw.ldfconst else object@ldfconst
    ret <- .mypredictLda(object@prior, levels(object@grp), ldf, ldfconst, x)
    if(ct)
        ret@ct <- mtxconfusion(object@grp, ret@classification)

    ret
})
