## Print function for ssanova objects
print.ssanova <- function(x,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    ## terms
    cat("Terms:\n")
    print.default(c(x$terms$labels,x$lab.p))
    cat("\n")
    ## terms overview
    cat("Number of unpenalized and penalized terms:\n\n")
    print.default(x$desc)
    cat("\n")
    if (x$method=="v") Method <- "GCV "
    if (x$method=="m") Method <- "GML.\n"
    if (x$method=="u") Method <- "Mallows CL "
    if (x$method=="m") cat("Smoothing parameters are selected by",Method)
    else cat("Smoothing parameters are selected by ",Method,"with alpha=",x$alpha,".",sep="")
    cat("\n")
    ## the rest are suppressed
    invisible()
}

## Print function for ssanova0 objects
print.ssanova0 <- function(x,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    ## terms
    cat("Terms:\n")
    print.default(c(x$terms$labels,x$lab.p))
    cat("\n")
    ## terms overview
    cat("Number of unpenalized and penalized terms:\n\n")
    print.default(x$desc)
    cat("\n")
    if (x$method=="v") Method <- "GCV.\n"
    if (x$method=="m") Method <- "GML.\n"
    if (x$method=="u") Method <- "Mallows CL.\n"
    cat("Smoothing parameters are selected by",Method)
    cat("\n")
    ## the rest are suppressed
    invisible()
}

## Print function for gssanova objects
print.gssanova <- function(x,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    ## terms
    cat("Terms:\n")
    print.default(c(x$terms$labels,x$lab.p))
    cat("\n")
    ## terms overview
    cat("Number of unpenalized and penalized terms:\n\n")
    print.default(x$desc)
    cat("\n")
    cat("Smoothing parameters are selected by CV with alpha=",x$alpha,".",sep="")
    cat("\n")
    ## the rest are suppressed
    invisible()
}

## Print function for ssden objects
print.ssden <- function(x,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    ## terms
    cat("Terms:\n")
    print.default(x$terms$labels)
    cat("\n")
    ## terms overview
    cat("Number of unpenalized and penalized terms:\n\n")
    print.default(x$desc)
    cat("\n")
    cat("Smoothing parameters are selected by CV with alpha=",x$alpha,".",sep="")
    cat("\n")
    ## the rest are suppressed
    invisible()
}

## Print function for sscden objects
print.sscden <- function(x,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    ## terms
    cat("Terms:\n")
    print.default(x$terms$labels)
    cat("\n")
    ## terms overview
    cat("Number of unpenalized and penalized terms:\n\n")
    print.default(x$desc)
    cat("\n")
    cat("Smoothing parameters are selected by CV with alpha=",x$alpha,".",sep="")
    cat("\n")
    ## the rest are suppressed
    invisible()
}

## Print function for sshzd objects
print.sshzd <- function(x,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    ## terms
    cat("Terms:\n")
    print.default(c(x$terms$labels,x$lab.p))
    cat("\n")
    ## terms overview
    cat("Number of unpenalized and penalized terms:\n\n")
    print.default(x$desc)
    cat("\n")
    cat("Smoothing parameters are selected by CV with alpha=",x$alpha,".",sep="")
    cat("\n")
    ## the rest are suppressed
    invisible()
}

## Print function for sshzd objects
print.sscox <- function(x,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    ## terms
    cat("Terms:\n")
    print.default(c(x$terms$labels,x$lab.p))
    cat("\n")
    ## terms overview
    cat("Number of unpenalized and penalized terms:\n\n")
    print.default(x$desc)
    cat("\n")
    cat("Smoothing parameters are selected by CV with alpha=",x$alpha,".",sep="")
    cat("\n")
    ## the rest are suppressed
    invisible()
}

## Print function for ssllrm objects
print.ssllrm <- function(x,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    ## terms
    cat("Terms:\n")
    print.default(x$terms$labels)
    cat("\n")
    ## terms overview
    cat("Number of unpenalized and penalized terms:\n\n")
    print.default(x$desc)
    cat("\n")
    cat("Smoothing parameters are selected by CV with alpha=",x$alpha,".",sep="")
    cat("\n")
    ## the rest are suppressed
    invisible()
}

## Print function for summary.ssanova objects
print.summary.ssanova <- function (x,digits=6,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n",sep="")
    cat("\nEstimate of error standard deviation:",x$sigma,"\n")
    ## residuals
    res <- x$res
    cat("\nResiduals:\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- structure(quantile(res), names = nam)
    print(rq,digits=digits)
    cat("Residual sum of squares:",x$rss)
    cat("\nR square:",x$r.squared)
    ## selected summaries
    cat("\n\nPenalty associated with the fit:",x$pen)
    cat("\n\n")
    invisible()
}

## Print function for summary.gssanova objects
print.summary.gssanova <- function (x,digits=6,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n",sep="")
    if (x$family%in%c("Gamma","inverse.gaussian")) {
        cat("\n(Dispersion parameter for ",x$family,
            " family estimated to be ",format(x$dispersion),")\n\n",sep="")
    }
    else {
        cat("\n(Dispersion parameter for ",x$family,
            " family taken to be ",format(x$dispersion),")\n\n",sep="")
    }
    ## residuals
    res <- x$res
    cat("Working residuals (weighted):\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- structure(quantile(res), names = nam)
    print(rq,digits=digits)
    cat("Residual sum of squares:",x$rss,"\n")
    ## deviance residuals
    res <- x$dev.res
    cat("\nDeviance residuals:\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- structure(quantile(res), names = nam)
    print(rq,digits=digits)
    cat("Deviance:",x$deviance)
    cat("\nNull deviance:",x$dev.null)
    ## selected summaries
    cat("\n\nPenalty associated with the fit:",x$pen)
    cat("\n\n")
    invisible()
}

## Print function for summary.gssanova objects
print.summary.gssanova0 <- function (x,digits=6,...)
{
    ## call
    cat("\nCall:\n",deparse(x$call),"\n",sep="")
    if (x$method=="u")
        cat("\n(Dispersion parameter for ",x$family,
            " family taken to be ",format(x$dispersion),")\n\n",sep="")
    if (x$method=="v")
        cat("\n(Dispersion parameter for ",x$family,
            " family estimated to be ",format(x$dispersion),")\n\n",sep="")
    ## residuals
    res <- x$res
    cat("Working residuals (weighted):\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- structure(quantile(res), names = nam)
    print(rq,digits=digits)
    cat("Residual sum of squares:",x$rss,"\n")
    ## deviance residuals
    res <- x$dev.res
    cat("\nDeviance residuals:\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- structure(quantile(res), names = nam)
    print(rq,digits=digits)
    cat("Deviance:",x$deviance)
    cat("\nNull deviance:",x$dev.null)
    ## selected summaries
    cat("\n\nPenalty associated with the fit:",x$pen)
    cat("\n\nNumber of performance-oriented iterations:",x$iter)
    cat("\n\n")
    invisible()
}
