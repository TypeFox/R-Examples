
##################################################################
### Function:  bde

##  bde: binned data density estimation

## histosmooth

## 2014/04/26: use regression method to estimate the parameters, also
## search the neiborhood of the initial estimates to find the MLEs.
## The method applies only to the top-coded data only.

## 2014/04/29:  add the other histosmooth algorithms to bde
## 2014/06/11: 

bde <-
    function(x, counts, nclass, breaks, bw,
             type="kde", from, to, gridsize=512L,
             lbound, conf.level)
    UseMethod("bde")

bde.default <-
    function(x, counts, nclass, breaks, bw,
             type="kde", from, to, gridsize=512L,
             lbound, conf.level)
{
    f.call <- match.call()
    xhist <- binning(x=x,counts=counts,breaks=breaks,bw=bw)
    out <- bde(xhist,type=type,from=from,to=to,
               gridsize=gridsize, lbound=lbound,
               conf.level=conf.level)
    out$call <- f.call
    out
}

bde.histogram <-
    function(x, counts, nclass, breaks, bw,
             type="kde", from, to, gridsize=512L,
             lbound, conf.level)
{
    xhist <- binning(x)
    out <- bde(xhist,type=type,from=from,to=to,
               gridsize=gridsize, lbound=lbound,
               conf.level=conf.level)
}

bde.bdata <-
    function(x, counts, nclass, breaks, bw,
             type="kde", from, to, gridsize=512L,
             lbound, conf.level)
{
    f.call <- match.call()
    ## support Dagum and Weibull only
    type <- match.arg(tolower(type),
                      c('kde', 'smkde','smoothkde',
                        'histospline','spline',
                        'lpr','npr','root-unroot',
                        'bootkde'))
    out <- switch(type,
                  'bootkde' = .bootkde(x,from=from,to=to,
                      gridsize=gridsize,conf.level=conf.level),
                  'lpr' = .histonpr(x,from=from,to=to,
                      gridsize=gridsize,conf.level=conf.level),
                  'npr' = .histonpr(x,from=from,to=to,
                      gridsize=gridsize,conf.level=conf.level),
                  'root-unroot' = .histonpr(x,from=from,to=to,
                      gridsize=gridsize, conf.level=conf.level),
                  'spline'=,'histospline'= .histospline(x, from=from,
                                to=to, gridsize=gridsize),
                  'smkde'=, 'smoothkde'=,
                  'kde' = .smkde(x,bandwidth=bw, from=from,to=to,
                      gridsize=gridsize)
                  )
    out$call <- f.call
    out
}

## out should be an R object "histosmooth": (1) y, (2) x, (3) lcb,
## (4) ucb/conf.level, (5) type, (6) xhist,  (8) pars (npar=.)


print.histosmooth <- function (x, digits = NULL, ...)
{
    cat("Call:  ", deparse(x$call), "\n", sep = "")
    print(summary(as.data.frame(x[c("x", "y")])), ...)
    cat("\n")
    invisible(x)
}


plot.histosmooth <- function (x, col=1, lwd=1, lty=1,
                         shade,border="gray",scb=FALSE,...)
{
    if(length(col)==1){
        col1 <- col; col2 <- col
    }else{
        col1 <- col[1]; col2 <- col[2]
    }
    if(length(lwd)==1){
        lwd1 <- lwd; lwd2 <- lwd
    }else{
        lwd1 <- lwd[1]; lwd2 <- lwd[2]
    }
    if(length(lty)==1){
        lty1 <- lty; lty2 <- lty
    }else{
        lty1 <- lty[1]; lty2 <- lty[2]
    }
    
    plot(x$x, x$y, col=col1, lty=lty1,lwd=lwd1,...)

    if(!is.null(x$ucb)&&!is.null(x$lcb)&&scb){
        if(missing(shade)){
            lines(x$ucb~x$x,col=col2,lty=lty2,lwd=lwd2,...)
            lines(x$lcb~x$x,col=col2,lty=lty2,lwd=lwd2,...)
        }else{
            y0 <- c(x$ucb, rev(x$lcb))
            x0 <- c(x$x, rev(x$x))
            polygon(x0, y0, col=shade, border=border,...)
            lines(x$x, x$y, col=col1,lty=lty1,lwd=lwd1,...)
        }
    }
    
    invisible(x)
}

lines.histosmooth <- function (x, col=1, lwd=1, lty=1,
                         shade,border="gray",scb=FALSE,...)
{
    if(length(col)==1){
        col1 <- col; col2 <- col
    }else{
        col1 <- col[1]; col2 <- col[2]
    }
    if(length(lwd)==1){
        lwd1 <- lwd; lwd2 <- lwd
    }else{
        lwd1 <- lwd[1]; lwd2 <- lwd[2]
    }
    if(length(lty)==1){
        lty1 <- lty; lty2 <- lty
    }else{
        lty1 <- lty[1]; lty2 <- lty[2]
    }
    
    if(!is.null(x$ucb)&&!is.null(x$lcb)&&scb){
        if(missing(shade)){
            lines(x$ucb~x$x,col=col2,lty=lty2,lwd=lwd2,...)
            lines(x$lcb~x$x,col=col2,lty=lty2,lwd=lwd2,...)
        }else{
            y0 <- c(x$ucb, rev(x$lcb))
            x0 <- c(x$x, rev(x$x))
            polygon(x0, y0, col=shade, border=border,...)
        }
    }
    lines(x$x, x$y, col=col1,lty=lty1,lwd=lwd1,...)
    
    invisible(x)
}



.histokde <-
    function(x,from=from,to=to,
             gridsize=gridsize, lbound=lbound)
    {
        0
    }
