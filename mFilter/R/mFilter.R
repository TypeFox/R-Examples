## Generic mFilter functions
## Part of mFilter package

mFilter <- function(x, ...) UseMethod("mFilter")

mFilter.default <- function(x, ...) mFilter.ts(x, ...)

mFilter.ts <- function(x, filter=c("HP","BK","CF","BW","TR"), ...)
{
    filt = match.arg(filter)
    call = match.call()
    ag = list(...)
    switch(filt,
           "HP" = {res <- hpfilter(x,freq=ag$freq,type=ag$type,drift=ag$drift)},
           "BK" = {res <- bkfilter(x,pl=ag$pl,pu=ag$pu,nfix=ag$nfix,type=ag$type,drift=ag$drift)},
           "CF" = {res <- cffilter(x,pl=ag$pl,pu=ag$pu,root=ag$root,drift=ag$drift,
                                 type=ag$type, nfix=ag$nfix,theta=ag$theta)},
           "BW" = {res <- bwfilter(x,freq=ag$freq,nfix=ag$nfix,drift=ag$drift)},
           "TR" = {res <- trfilter(x,pl=ag$pl,pu=ag$pu,drift=ag$drift)}
                  )
    res$xname <- deparse(substitute(x))
    return(res)
}

print.mFilter <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    if (!inherits(x, "mFilter"))
    stop("method is only for mFilter objects")

    # Title:
    cat("\nTitle:\n ")
    cat(x$title, "\n")

    ## Call:
    cat("\nCall:\n ")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

    ## Method:
    cat("\nMethod:\n ", x$method, "\n", sep = "")

    ## Filter Type:
    cat("\nFilter Type:\n ", x$type, "\n", sep = "")

    ## Series
    cat("\nSeries:\n ", x$xname, "\n\n", sep = "")

    names <- c(x$xname,"Trend","Cycle")
    out <- cbind(x$x,x$trend,x$cycle)
    colnames(out) <- names
    rownames(out) <- time(x$x)
    if(any(frequency(x$x) == c(4,12)))
        print(out, digits=digits)
    else
        print(as.data.frame(out),digits = digits)

    invisible(x)
}

summary.mFilter <- function(object, digits = max(3, getOption("digits") - 3), ...)
{
    if (!inherits(object, "mFilter"))
        stop("method is only for mFilter objects")

    ## Title:
    cat("\nTitle:\n ")
    cat(object$title, "\n")

    ## Call:
    cat("\nCall:\n ")
    cat(paste(deparse(object$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

    ## Method:
    cat("\nMethod:\n ", object$method, "\n", sep = "")

    ## Filter Type:
    cat("\nFilter Type:\n ", object$type, "\n", sep = "")

    ## Series
    cat("\nSeries:\n ", object$xname, "\n", sep = "")

    names <- c(object$xname,"Trend","Cycle")
    out <- cbind(object$x,object$trend,object$cycle)
    colnames(out) <- names
    rownames(out) <- time(object$x)

    cat("\nDescriptive Statistics:\n ", "\n", sep = "")
    print(summary(out), digits = digits)
    #browser()
    gof <- function(object)
    {
        res <- object$cycle
        pe <- res/object$x
        out <- c(mean(res,na.rm=TRUE), mean(res^2,na.rm=TRUE),
                 mean(abs(res),na.rm=TRUE), mean(pe,na.rm=TRUE), mean(abs(pe),na.rm=TRUE))
        names(out) <- c("ME","MSE","MAE","MPE","MAPE")
        return(out)
    }
    cat("\nIn-sample error measures:\n")
    print(gof(object), digits = digits)

    cat("\n")
    #if(any(frequency(object$x) == c(4,12)))
        #print(out)
    #else
        #print(as.data.frame(out))

    invisible(object)
}


plot.mFilter <- function(x, reference.grid = TRUE, col = "steelblue", ask=interactive(), ...)
{
    if (!inherits(x, "mFilter"))
        stop("method is only for mFilter objects")
    opar <- par(no.readonly=TRUE)
    par(ask=ask,mfrow=c(2,1),mar=c(3,2,2,1))

    ag <- list(...)
    if(is.null(ag$main))
        main <- paste(x$title, "of", x$xname)

    ylim <- range(c(x$x,x$trend),na.rm=TRUE)
    plot(x$x,type="l",main=main,ylab="",ylim=ylim,col=col,...)
    lines(x$trend,col="red")
    if (reference.grid) grid()
    legend("topleft",legend=c(x$xname,"trend"),col=c(col,"red"),lty=rep(1,2),ncol=2)
    plot(x$cycle,type="l",main="Cyclical component (deviations from trend)",
         ylab="",col=col,...)
    if (reference.grid) grid()
    par(opar)

    invisible(x)
}

residuals.mFilter <- function(object, ...) return(object$cycle)

fitted.mFilter <- function(object, ...) return(object$trend)




