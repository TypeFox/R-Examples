### Methods for Weibull or Weibull-related distributions

## Created on 2014/06/12 


##################################################
##  Weibull(x,kappa, lambda)
##  f(x) = k/l*(x/l)^(k-1)*exp(-(x/l)^k)
##  in R Weibill(x,kappa,lambda)
##       -- kappa==shape, lambda=scale
##################################################

fit.Weibull <- function(x, dist="Weibull")
    UseMethod("fit.Weibull")

fit.Weibull.default <- function(x, dist="Weibull")
{
    out0 <- .WeibullPlot(x)#initial values
    xhist <- hist(x,plot=FALSE)
    dist <- match.arg(tolower(dist),
                      c("weibull","gwd","ewd"))
    res <- switch(dist,
                  'weibull' = .weibullMLE(x,out0),
                  'gwd' = .gwdBMLE(x),
                  'ewd' = .ewdBMLE(x,out0))
    structure(
        list(xhist=xhist, pars=res,dist=dist,
             par0=out0, call = match.call()),
        class="GWD")
}

fit.Weibull.bdata <- function(x,dist='Weibull')
{
    out0 <- .WeibullbPlot(x)#initial values
    dist <- match.arg(tolower(dist),
                      c("weibull","gwd","ewd"))
    res <- switch(dist,
                  'weibull' = .weibullBMLE(x,out0),
                  'gwd' = .gwdBMLE(x),
                  'ewd' = .ewdBMLE(x,out0))
    structure(
        list(xhist=x, pars=res,dist=dist,
             par0=out0, call = match.call()),
        class="GWD")
}

fit.Weibull.histogram <- function(x,dist='Weibull')
{
    xhist <- binning(x)
    fit.Weibull(xhist,dist=dist)
}

print.GWD <- function(x,...){
    if(x$dist == "weibull"){
        cat("The Weibull Distribution\n  Call:  ",
            deparse(x$call), sep = "")
        if(is.null(x$pars)){
            warning("Weibull fitting failed!")
        }else{
            a <- x$pars[1]
            b <- x$pars[2]
            cat("\n  Parameters: = (", a, ", ",
                b, ")\n", sep='')
            cat("  Range: =[ 0, Inf)\tMedian: =",
                qweibull(0.5, a,b), "\n",sep='')
        }
    }else if(x$dist=='ewd'){
        cat("The Exponentiated Weibull Distribution\n  Call:  ",
            deparse(x$call), sep = "")
        if(is.null(x$pars)){
            warning("EWD fitting failed!")
        }else{
            a <- x$pars[1]
            b <- x$pars[2]
            p <- x$pars[3]
            cat("\n  Parameters: = (", a, ", ",
                b, ", ", p,")\n", sep='')
            cat("  Range: =[ 0, Inf)\tMedian: =",
                qewd(0.5, a,b,p), "\n",sep='')
        }
    }else if(x$dist=='gwd'){
        cat("The Generalized Weibull Distribution\n  Call:  ",
            deparse(x$call), sep = "")
        if(is.null(x$pars)){
            warning("GWD fitting failed!")
        }else{
            a <- x$pars[1]
            b <- x$pars[2]
            p <- x$pars[3]
            cat("\n  Parameters: = (", a, ", ",
                b, ", ", p,")\n", sep='')
            cat("  Range: =[ 0, Inf)\tMedian: =",
                qgwd(0.5, a,b,p), "\n",sep='')
        }
    }
}

plot.GWD <- function(x,hist=TRUE,...){
    if(is.null(x$pars))
        warning("GLD fitting failed!")
    else{
        res <- .bdataTohist(x$xhist)
        k <- length(res$mids)
        a <- res$breaks[1]
        b <- res$breaks[k+1]
        x0 <- seq(a, b, length=100)

        if(x$dist == "weibull"){
            a <- x$pars[1]
            b <- x$pars[2]
            f0 <- dweibull(x0, a,b)
        }else if(x$dist=='ewd'){
            a <- x$pars[1]
            b <- x$pars[2]
            p <- x$pars[3]
            f0 <- dewd(x0, a,b,p)
        }else if(x$dist=='gwd'){
            a <- x$pars[1]
            b <- x$pars[2]
            p <- x$pars[3]
            f0 <- dgwd(x0, a,b,p)
        }

        if(hist){
            plot(res,...)
            lines(f0~x0,...)
        }else{
            plot(f0~x0, ...)
        }
    }
}

lines.GWD <- function(x,...){
    if(is.null(x$pars))
        warning("GLD fitting failed!")
    else{
        res <- .bdataTohist(x$xhist)
        k <- length(res$mids)
        a <- res$breaks[1]
        b <- res$breaks[k+1]
        x0 <- seq(a, b, length=100)

        if(x$dist == "weibull"){
            a <- x$pars[1]
            b <- x$pars[2]
            f0 <- dweibull(x0, a,b)
        }else if(x$dist=='ewd'){
            a <- x$pars[1]
            b <- x$pars[2]
            p <- x$pars[3]
            f0 <- dewd(x0, a,b,p)
        }else if(x$dist=='gwd'){
            a <- x$pars[1]
            b <- x$pars[2]
            p <- x$pars[3]
            f0 <- dgwd(x0, a,b,p)
        }

        lines(f0~x0,...)
    }
}


.WeibullbPlot <- function(x){
    xhist <- x
    x <- xhist$breaks[-1]
    n <- length(x)
    N <- sum(xhist$counts)
    Fhat <- cumsum(xhist$counts)/N
    Fhat[n] <- (N-0.3)/(N+0.4)
    ##    Fhat <- (rank(x)-0.3)/(n+0.4)
    
    ly <- log(-log(1-Fhat))

    tol <- 1e-10
    if(xhist$top.coded==1||xhist$top.coded==2)
        x[n] <- 1e8
    lx <- log(x+tol)
    out <- lm(ly~lx)
    kappa <- out$coef[[2]]
    lambda <- exp(-out$coef[[1]]/kappa)
    c(kappa, lambda)
}

.WeibullPlot <- function(x){
    ## x must be raw data, not histogram or bdata or others
    stopifnot(class(x)=="numeric")
    n <- length(x)
    Fhat <- (rank(x)-0.3)/(n+0.4)
    ly <- log(-log(1-Fhat))

    tol <- 1e-10
    x.inf <- !is.finite(x)
    if(any(x.inf)){
        x[x.inf] <- sign(x[x.inf])/tol
    }
    lx <- log(x+tol)
    out <- lm(ly~lx)
    kappa <- out$coef[[2]]
    lambda <- exp(-out$coef[[1]]/kappa)
    c(kappa, lambda)
}

.weibull.k <- function(k,x){
    x.out <- .check(x,sign='nn')
    x <- x.out$x
    stopifnot(k>0)
    xmin <- min(x)
    lxmin <- log(xmin)
    sum1 <- sum(x^k*log(x))
    sum2 <- sum(x^k)
    sum3 <- mean(log(x))
    sum1/sum2-sum3-1/k
}

.weibullBMLE <- function(x, pars)
{
    tol <- 0.0000001

    .weibull.bllk <- function(pars){
        res <- -sum(counts * log(diff(pweibull(breaks, pars[1], pars[2]))))
        if(!is.finite(res)) res <- 99999999999.
        res
    }

    breaks <- x$breaks
    counts <- x$counts

    out <- optim(pars, .weibull.bllk, ## gr=.weibull.grbllk,
                 method="L-BFGS-B", lower=c(tol,tol))
    as.numeric(out$par)
}


.weibullMLE <- function(x, pars)
{
    kappa <- uniroot(.weibull.k,c(pars[1]*.80,pars[1]*1.80),x=x)$root
    lambda <- (mean(x^kappa))^(1/kappa)
    c(kappa,lambda)
}


.weibull.grbllk <- function(pars,breaks, counts){
    tmp1 <- diff(pweibull(breaks, pars[1], pars[2]))
    tmp2 <- dweibull(breaks, pars[1], pars[2])
    tmp2k <- .weibull.grk(breaks, pars[1], pars[2])
    tmp2l <- .weibull.grl(breaks, pars[1], pars[2])
    dk <- sum(counts * diff(tmp2*tmp2k)/tmp1)
    dl <- sum(counts * diff(tmp2*tmp2l)/tmp1)
    c(dk,dl)
}

.weibull.grk <- function(x,kappa,lambda){
    f <- dweibull(x,kappa,lambda)
    tmp <- f/kappa+log(x)-log(kappa)-(x/lambda)^kappa*(log(x)-log(kappa))
    tmp * f
}

.weibull.grl <- function(x,kappa,lambda){
    f <- dweibull(x,kappa,lambda)
    tmp <- ((x/lambda)^kappa-1)*kappa/lambda
    tmp * f
}

### EWD = Weibull^alpha:  consider only for grouped data

.ewdBMLE <- function(x, pars)
{
    if(class(x) != 'bdata'){
        warning("Data 'x' was grouped")
        xhist <- hist(x, plot=FALSE)
        x <- binning(xhist)
    }

    .ewd.bllk <- function(pars){
        res <- -sum(counts * log(diff(pewd(breaks, pars[1], pars[2], pars[3]))))
        if(!is.finite(res)) res <- 99999999999.
        res
    }

    breaks <- x$breaks
    counts <- x$counts
    
    tol <- 0.0000001
    pars <- c(1,pars)
    out <- optim(pars, .ewd.bllk, ## gr=.weibull.grbllk,
                 method="L-BFGS-B", lower=c(tol,tol, tol))
    as.numeric(out$par)
}


## GWD: for grouped data only

.gwdBMLE <- function(x)
{
    if(class(x) != 'bdata'){
        warning("Data 'x' was grouped")
        xhist <- hist(x, plot=FALSE)
        x <- binning(xhist)
    }

    .gwd.bllk <- function(pars){
        res <- -sum(counts * log(diff(pgwd(breaks, pars[1], pars[2], pars[3]))))
        if(!is.finite(res)) res <- 99999999999.
        res
    }

    breaks <- x$breaks
    counts <- x$counts

    tol <- 0.0000001
    x.sum <- .summary.bdata(x)
    ## try the exact quantiles first
    sigma <- sqrt(x.sum$var)
    qlevels <- x.sum$exact$levels
    qtls <- x.sum$exact$qtls
    pars <- .gwd.pm(sigma, qlevels, qtls)
    if(is.null(pars)){
        sigma <- sqrt(x.sum$var)
        qlevels <- x.sum$approx$levels
        qtls <- x.sum$approx$qtls
        pars <- .gwd.pm(sigma, qlevels, qtls)
    }
    if(is.null(pars)) pars <- c(1,1)
    out <- optim(pars, .gwd.bllk, ## gr=.weibull.grbllk,
                 method="L-BFGS-B", lower=c(tol,tol,-Inf))
    as.numeric(out$par)
}

.gwd.pm <- function(sigma, qlevels, qtls){
    ## to estimate the parameters of GWD using percentile matching
    ## method.

    .fgwdpm <- function(lambda){
        abs(log((1-(1-p2)^lambda)/lambda)/log((1-(1-p1)^lambda)/lambda)-b)
    }


    x.na <- is.na(qtls)
    k <- sum(!x.na)
    if(k > 1){
        qtls <- qtls[!x.na]
        qlevels <- qlevels[!x.na]
        if(k >= 4){
            q1 <- qtls[2]; q2 <- qtls[4]
            p1 <- qlevels[2]; p2 <- qlevels[4];
        }else if(k == 3){
            q1 <- qtls[1]; q2 <- qtls[3]
            p1 <- qlevels[1]; p2 <- qlevels[3];
        }else{
            q1 <- qtls[1]; q2 <- qtls[2]
            p1 <- qlevels[1]; p2 <- qlevels[2];
        }
        
        b <- (log(q2)-log(sigma))/(log(q1)-log(sigma))
        pars <- -0.5
        out <- optim(pars, .fgwdpm, method="L-BFGS-B")
        lambda <- as.numeric(out$par)
        if(lambda != 0)
            alpha <- log(q2/sigma)/log((1-(1-p2)^lambda)/lambda)
        else
            alpha <- -log(q2/sigma)/log(1-p2)
        res <- c(alpha, sigma, lambda)

    }else{
        res <- NULL
    }
    res
}


