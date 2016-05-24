## Created on 06/18/2014.

gof <- function(x, dist, pars, method='ks')
    UseMethod("gof")

gof.default <- function(x, dist, pars, method='ks')
{
    if(!is.character(method)) method <- 'ks'
    
    method <- match.arg(tolower(method),
                        c("chisquare","ks","kolmogorov-smirnov"))
    
    if(method=='chisquare'){
        xhist <- hist(x, plot=FALSE)
        out <- .gof.chi(x,dist=dist, pars=pars)
    }else{
        out <- .gof.ks(x,dist=dist, pars=pars)
    }
    out
}

gof.bdata <- function(x, dist, pars, method='ks')
{
    if(!is.character(method)) method <- 'ks'
    
    method <- match.arg(tolower(method),
                        c("chisquare","ks","kolmogorov-smirnov"))

    if(method=='chisquare'){
        out <- .gof.chi(x,dist=dist, pars=pars)
    }else{
        out <- .gof.ks(x,dist=dist, pars=pars)
    }
    out
}

.gof.chi <- function(x, dist, pars){

    if(missing(dist)) dist <- 'normal'
    dist <- match.arg(tolower(dist),## add RST-GLD later
                      c("normal","weibull","fmkl","gb"))
    out <- .compact(x)
    cts <- out$counts; xn <- out$breaks
    
    Fx <- switch(dist,
                 "normal" = pnorm(xn, pars),
                 "weibull" = pweibull(xn, pars),
                 "fmkl" = pgld(xn, pars),
                 "gb" = pgbeta(xn, pars)
                 )
    Nhat <- sum(cts) * diff(Fx)
    X2 <- sum((Nhat - cts)^2/Nhat)
    ps <- x$nclass - length(pars) - 1
    parameters <- ifelse(ps>0, ps, 1)
    pval <- pchisq(X2, parameters, lower.tail=FALSE)
    dname <- x$name
    
    names(X2) <- "Chi-square statistic"
    names(parameters) <- "df"
    RVAL <- list(statistic = X2, parameter = parameters,
                 method = "Chi-square Goodness-of-Fit Test",
                 alternative = "two-sided",
                 p.value= pval, observed = cts,
                 expected = Nhat, residuals = sqrt(X2),
                 data.name=dname)
    class(RVAL) <- "htest"
    return(RVAL)
}


.gof.ks <- function(x, dist, pars){

    if(class(x)=="bdata"){
        xn <- x$breaks
        Fn <- c(0,x$counts/sum(x$counts))
    }else{
        x.t <- table(x)
        xn <- as.numeric(x.t)
        Fn <- x.t/sum(x.t)
        Fn <- c(0, as.numeric(Fn))
    }

    if(missing(dist)) dist <- 'normal'
    dist <- match.arg(tolower(dist),## add RST-GLD later
                      c("normal","weibull","fmkl","gb"))
    Fx <- switch(dist,
                 "normal" = pnorm(xn, pars),
                 "weibull" = pweibull(xn, pars),
                 "fmkl" = pgld(xn, pars),
                 "gb" = pgbeta(xn, pars)
                 )

    D <- max(abs(Fx - Fn))
    pval <- .Fortran(.F_KSPvalue, pv=as.double(D))$pv
    dname <- x$name
    
    names(D) <- "D"
    RVAL <- list(statistic = D,
                 method = "One-sample Kolmogorov-Smirnov Test",
                 alternative = "two-sided",
                 p.value= pval,
                 data.name=dname)
    class(RVAL) <- "htest"
    return(RVAL)
}

gof.FMKL <- function(x, dist, pars, method='ks')
{
    xhist <- x$xhist
    gof(xhist, dist="fmkl", pars=x$pars, method=method)
}

gof.histogram <- function(x, dist, pars, method='ks')
{
    xhist <- binning(x)
    gof(xhist, dist=dist, pars=pars, method=method)
}

gof.GB <- function(x, dist, pars, method='ks')
{
    xhist <- x$xhist
    gof(xhist, dist="gb", pars=x$pars, method=method)
}
