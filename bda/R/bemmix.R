## 2014-04-10: add an option that the number of components are unknown

## x of size nr+1, sorted, class limits; y: frequencies/counts, of
## size nr.  mu: initial mean values of the components.

## k is the number of parameters; L is the likelihood; n is the sample
## size

## AIC = 2*k - 2*ln(L)
## BIC = k*[ln(n)+ln(2*pi)] - 2*ln(L)
## BIC = k*ln(n) - 2*ln(L) => for large n
## AICc = AIC + 2*k*(k+1)/(n-k-1)

## 2014/04/19: if 'k' is not specified, we try k=1,2,...,10 and choose
## the best fit using AICc.

fit.mixnorm <-
    function(x,k,mu,s,p, x.range, lognormal=FALSE)
    UseMethod("fit.mixnorm")

fit.mixnorm.default <-
    function(x,k,mu,s,p, x.range,lognormal=FALSE){
        f.call <- match.call()
        xhist <- hist(x, plot=FALSE)
        res <- fit.mixnorm(xhist,k=k,mu=mu,s=s,p=p,
                        x.range=x.range,lognormal=lognormal)
        res$call <- f.call
}

fit.mixnorm.histogram <-
    function(x, k, mu, s, p, x.range, lognormal=FALSE){
        f.call <- match.call()
        xhist <- binning(x)
    
        res <- fit.mixnorm(x=xhist,
                        k=k, mu=mu,s=s,p=p,
                        x.range=x.range,
                        lognormal=lognormal)
        res$call <- f.call
        res
    }

fit.mixnorm.bdata <-
    function(x, k, mu, s, p, x.range, lognormal=FALSE){
        f.call <- match.call()

        res <- .cemmix(breaks=x$breaks,freq=x$counts,
                       x$top.coded, mu=mu,s=s,p=p,k=k,
                       x.range=x.range, lognormal=lognormal)
        res$xhist <- x
        res$call <- f.call

        xpts <- c(res$x.range[1],x$breaks[1],
                  x$breaks[x$nclass+1],res$x.range[2])
        Fx <- cdf.nmix(res, x0=xpts)
        res$nl <- (Fx$y[2]-Fx$y[1])*sum(x$counts)
        res$nu <- (Fx$y[4]-Fx$y[3])*sum(x$counts)
        res
    }


.cemmix <- function(breaks,freq,top.coded, mu,s,p,k,
                    x.range,lognormal=FALSE){

    stopifnot(is.numeric(breaks))
    stopifnot(is.numeric(freq))
    stopifnot(!any(is.na(breaks)))
    stopifnot(!any(is.na(freq)))

    freq <- round(freq) # in case not integer(s)
    if(any(freq < 0)) stop("Counts cannot be negative")

    ## breaks should be monotone increasing
    stopifnot(any(diff(breaks) > 0))

    nr <- length(freq)
    stopifnot(nr>3)
    stopifnot(length(breaks) == nr +1)

    ## McLachlan & Jones (1988): in the case of not truncated
    ## data, the equations require two extra classes (-\infty, x0)
    ## and (xr,+\infty) with corresponding frequencies n_{l} and
    ## n_{u}.
    itrunc <- 0; trunc <- FALSE # output
    nl <- 0; nu <- 0;

    if(missing(x.range)) x.range <- c(-Inf, Inf)

    if(top.coded == 2){
        itrunc <- 1;
        trunc <- TRUE;
        nl <- freq[1];
        nu <- freq[nr];
        breaks <- breaks[-c(1,nr+1)];
        freq <- freq[-c(1,nr)];
        x.range <- c(-Inf,Inf)
    }

    if(top.coded == 1){
        itrunc <- 1;
        trunc <- TRUE;
        nl <- 0
        nu <- freq[nr];
        breaks <- breaks[-(nr+1)];
        freq <- freq[-nr];
        x.range[2] <- Inf
    }

    if(top.coded == -1){
        itrunc <- 1;
        trunc <- TRUE;
        nl <- freq[1];
        nu <- 0
        breaks <- breaks[-1];
        freq <- freq[-1];
        x.range[1] <- -Inf
    }
    
    ## renew the number of groups (bins)
    nr <- length(freq)

    ## whether to log-transform data
    c.glog <- 0;  # a constant for glog-transformation
    if(lognormal){
        c.glog <- ifelse(min(breaks)<0.25, 0.25-min(breaks), min(breaks))
        breaks <- log(breaks + c.glog)
    }
    
    ## Apply EM for model fitting
    x0 <- breaks[1]; x1 <- breaks[-1]
    xc <- (breaks[-1]+breaks[-(nr+1)])*0.5

    ## if both 'mu' and 'k' are given, first use 'mu'
    if(missing(mu)){
        if(missing(k))
            ng <- 3
        else
            ng <- round(k)
        ng <- min(ng, nr)
        mu <- sample(xc, size=ng, prob = freq/sum(freq))
    }else{
        if(missing(k))
            ng <- length(mu)
        else
            stopifnot(length(mu)==k)
    }
    
    if(missing(p)) p <- rep(1/ng,ng)  # initial proportions
    else{
        stopifnot(!any(p<=0))
        p <- p/sum(p)
    }
    
    if(missing(s)) v <- rep(var(rep(xc, freq)), ng)
    else v <- s^2
    
    stopifnot(length(p) == ng)
    stopifnot(length(v) == ng)

    ## input as WK (workspace size) and output as error info
    wk <- max(500, 3*ng*(nr+2+3*ng)+(nr+2))

    out <- .Fortran(.F_emmix,
                    as.integer(freq), as.integer(nr), as.integer(ng),
                    as.double(x0), as.double(x1),
                    p=as.double(p), mu=as.double(mu), v=as.double(v),
                    xlogl=double(1), ifault=as.integer(wk),
                    iter = as.integer(itrunc), 
                    nl=as.integer(nl), nu=as.integer(nu))
    llk0 <- out$xlogl
    K <- ifelse(ng == 1, 2., 3. * ng - 1.);
    N <- sum(freq);
    xspan <- diff(range(breaks)) 
    AIC = -2.0 * llk0 + 2.0 * K;
    AICc = AIC + 2.0 * K * (K+1.) / (N - K - 1.0);
    BIC = -2.0 * llk0 + log(N*2.0*pi) * K;

    structure(list(xhist=NULL, 
                   ng=ng,p=out$p, mu=out$mu-c.glog, s=sqrt(out$v),
                   llk = out$xlogl, nl = nl, nu = nu,
                   iter = out$iter, ifault = out$ifault,
                   "AIC"=AIC, "BIC"=BIC,"AICc"=AICc,
                   lognormal = lognormal, c.glog=c.glog,
                   x.range = x.range,
                   call = match.call(),
                   trunc = trunc),
              class='nmix')
}

print.nmix <- function(x,...)
    {
        cat("\nCall:\n\t", deparse(x$call), "\n")
        tmp.names <- paste("Comp.",c(1:length(x$mu)),sep='')
        tmp <- cbind(x$p, x$mu, x$s)
        row.names(tmp) <- tmp.names
        tmp <- t(tmp)
        row.names(tmp) <- c("Proportion","Mean", "Std.Dev")
        cat("\nParameter estimates:\n")
        print(tmp, ...)
        if(!x$trunc){
            cat("\nEstimates of the lower and upper class frequencies:\n")
            cat("Lower class: nl=", x$nl, "\tUpper class: nu=", x$nu)
        }
        cat("\n")
    }

plot.nmix <- function(x,hist=TRUE,...){
    res <- .bdataTohist(x$xhist)
    out <- pdf.nmix(x)
    
    if(hist){
        plot(res,...)
        lines(out,...)
    }else{
        plot(out, ...)
    }
}

lines.nmix <- function(x,...){
    out <- pdf.nmix(x)
    lines(out,...)
}

pdf.nmix <- function(x, x0, from, to, gridsize=512){
    stopifnot(class(x)=='nmix')
    if(missing(x0)){
        xhist <- x$xhist
        xhist2 <- .bdataTohist(xhist)
        k <- xhist$nclass
        if(missing(from))
            from <- xhist2$breaks[1]
        if(missing(to))
            to <- xhist2$breaks[k+1]
        stopifnot(to > from)
        x0 <- seq(from, to, length=gridsize)
    }
    x.range <- x$x.range
    if(x$lognormal){
        x0 <- log(x0 - min(x0))
        x.range <- log(x.range)
    }
    y0 <- dmixnorm(x0,p=x$p, mean=x$mu, sd=x$s)
    tot.mass <- diff(pmixnorm(x.range,p=x$p, mean=x$mu, sd=x$s))
    y0 <- y0/tot.mass
    list(x = x0, y = y0);
}

cdf.nmix <- function(x, x0, from, to, gridsize=512){
    stopifnot(class(x)=='nmix')
    if(missing(x0)){
        xhist <- x$xhist
        xhist2 <- .bdataTohist(xhist)
        k <- xhist$nclass
        if(missing(from))
            from <- xhist2$breaks[1]
        if(missing(to))
            to <- xhist2$breaks[k+1]
        stopifnot(to > from)
        x0 <- seq(from, to, length=gridsize)
    }
    x.range <- x$x.range
    if(x$lognormal){
        x0 <- log(x0 - min(x0))
        x.range <- log(x.range)
    }
    y0 <- pmixnorm(x0,p=x$p, mean=x$mu, sd=x$s)
    tot.mass <- diff(pmixnorm(x.range,p=x$p, mean=x$mu, sd=x$s))
    y0 <- y0/tot.mass
    list(x = x0, y = y0);
}

.dnmix <- function(res,x, x0){
    ##  this subroutine is not used.  Use pdf.nmix instead
    if(missing(x0)){
        k <- length(res$mids)
        a <- res$breaks[1]
        b <- res$breaks[k+1]
        gridsize <- 100
        x0 <- seq(a, b, length=gridsize)
    }
    x.range <- x$x.range
    if(x$lognormal){
        x0 <- log(x0 - min(x0))
        x.range <- log(x.range)
    }
    y0 <- dmixnorm(x0,p=x$p, mean=x$mu, sd=x$s)
    tot.mass <- diff(pmixnorm(x.range,p=x$p, mean=x$mu, sd=x$s))
    y0 <- y0/tot.mass
    list(x = x0, y = y0);
}


