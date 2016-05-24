###  To fit a smooth curve to histogram-type data using GB
###  (4-parameter Beta distribution)

## Created on 2014/06/18

fit.GB <- function(x, lbound, ubound)
{
    if(class(x) == 'histogram')
        x <- binning(x)
    if(class(x) != 'bdata')
        stop("The current version support histogram-type data only.")

    if(x$top.coded !=0)
        stop("Unbounded support")
    
    if(missing(lbound)){
        a <- x$breaks[1]
    }else if(is.numeric(lbound)){
        if(lbound > x$breaks[1])
            stop("'lbound' larger than the lower limit of the first class")
        if(is.finite(lbound))
            a <- lbound
        else
            stop("invalid 'lbound' value")
    }else stop("'lbound' is not numeric")

    if(missing(ubound)){
        b <- x$breaks[x$nclass+1]
    }else if(is.numeric(ubound)){
        if(ubound < x$breaks[x$nclass+1])
            stop("'ubound' smaller than the upper limit of the last class")
        if(is.finite(ubound))
            b <- ubound
        else
            stop("invalid 'ubound' value")
    }else stop("'ubound' is not numeric")

    ## compute initial values
    breaks <- x$breaks
    counts <- x$counts
    k <- x$nclass
    mids <- 0.5*(breaks[-1]+breaks[-k])
    mids <- (mids-a)/(b-a)
    xbar <- sum(mids * counts)/sum(counts)
    xvar <- sum((mids-xbar)^2*counts)/(sum(counts)-1)

    alpha <- (xbar*(1-xbar)/xvar -1)*xbar
    lambda <- (1/xbar - 1)*alpha
    
    res <- .gb.mle(x, alpha, lambda)
    pars <- c(a,b,res)
    
    structure(
        list(xhist=x, pars=pars,method="MLE",
             call = match.call()),
        class="GB")
}

print.GB <- function(x,...){
    cat("GB\n  Call:  ", deparse(x$call), sep = "")
    if(is.null(x$pars)){
        warning("GB fitting failed!")
    }else{
        pars <- x$pars
        cd <- pars[3] + pars[4]
        ba <- pars[2] - pars[1]
        cat("\n  Parameters: = (", pars[1], ", ",
            pars[2], ", ", pars[3], ", ",
            pars[4],")\n", sep='')
        cat("  Range: =[",pars[1],", ", pars[2],"]\tMean: =",
            pars[3]/cd*ba+pars[1],
            "\tVar: =",
            ba^2*pars[3]*pars[4]/(cd^2*(cd+1)),
            "\n",sep='')
    }
}

plot.GB <- function(x,hist=TRUE,...){
    if(is.null(x$pars))
        warning("GLD fitting failed!")
    else{
        res <- .bdataTohist(x$xhist)
        k <- length(res$mids)
        a <- res$breaks[1]
        b <- res$breaks[k+1]
        x0 <- seq(a, b, length=100)
        y0 <- (x0-x$pars[1])/(x$pars[2]-x$pars[1])
        f0 <- dbeta(y0, x$pars[3], x$pars[4])/(x$pars[2]-x$pars[1])
        if(hist){
            plot(res,...)
            lines(f0~x0,...)
        }else{
            plot(f0~x0, ...)
        }
    }
}

lines.GB <- function(x,...){
    if(is.null(x$pars))
        warning("GB fitting failed!")
    else{
        res <- .bdataTohist(x$xhist)
        k <- length(res$mids)
        a <- res$breaks[1]
        b <- res$breaks[k+1]
        x0 <- seq(a, b, length=100)
        y0 <- (x0-x$pars[1])/(x$pars[2]-x$pars[1])
        f0 <- dbeta(y0, x$pars[3], x$pars[4])/(x$pars[2]-x$pars[1])
        lines(f0~x0,...)
    }
}

    

.gb.mle <- function(x,alpha,lambda){
    .g.gbmle <- function(lmd){
        breaks <- x$breaks
        counts <- x$counts
        Fx <- diff(pbeta(breaks, lmd[1], lmd[2]))
        sele <- Fx <= 0
        Fx[sele] <- 1e-10
        -sum(counts * log(Fx))
    }


    ## find the MLEs
    if(alpha <= 0 || lambda <= 0){
        pars <- c(1,1)
    }else{
        pars <- c(alpha, lambda)
    }
    
    out <- optim(pars, .g.gbmle, method="L-BFGS-B",
                 lower = c(0, 0), upper = c(Inf,Inf))
    pars <- as.numeric(out$par)
}

