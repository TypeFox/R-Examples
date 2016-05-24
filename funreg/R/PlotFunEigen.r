#' @title plot method for funeigen object
#' @description Creates a visual representation of some of the 
#' information in an object of class \code{funeigen} 
#' (i.e., in an eigenfunction decomposition of a 
#' functional variable).  Several kinds of plots are 
#' available.
#' @param x An object of class \code{funeigen} 
#' @param type A character string telling the kind 
#' of information to include in the plot.  It may be \code{functions},
#' \code{eigenfunctions}, \code{eigenvalues}, \code{mean}, \code{covariance}, or \code{correlation}.
#' @param how.many How many fitted curves to show (in a plot of 
#' fitted curves; the default is all of them), or how many estimated 
#' eigenfunctions to show (in a plot of eigenfunctions; the default is
#' all of them)
#' @param xlab Label for the x axis of the plot.
#' @param ylab Label for the y axis of the plot.
#' @param ... Other optional arguments to be passed on to the plot function.
#' @export
#' @S3method plot funeigen
#' @method plot funeigen
plot.funeigen <- function(x,
                          type="correlation",
                          how.many=NULL, 
                          xlab="",
                          ylab="",
                          ...) {
    type <- tolower(type);
    stopifnot(class(x)=="funeigen");
    fit <- fitted(x,type=type);
    mids <- fitted(x,type="midpoints");
    if (type=="functions" | type=="centered") {
        if (is.null(how.many)) {how.many <- nrow(fit);}
        plot(x=mids,
             y=0*mids,
             type="n",
             xlim=c(min(mids),max(mids)),
             ylim=c(min(fit),max(fit)),
             xlab=xlab,ylab=ylab, ...);
        for (i in 1:nrow(fit)) {
            lines(mids,fit[i,]);
            if (how.many < 10) {text(mids,fit[i,],labels=i);}
        }
    }
    if (type=="eigenfunctions") {
        if (is.null(how.many)) {how.many <- ncol(fit);}
        plot(x=mids,
             y=0*mids,
             type="n",
             xlim=c(min(mids),max(mids)),
             ylim=c(min(fit),max(fit)),
             xlab=xlab,ylab=ylab, ...);
        for (i in 1:how.many) {
            lines(mids,fit[,i]);
            if (how.many < 10) {text(mids,fit[,i],labels=i);}
        }
    }
    if (type=="eigenvalues") {
        plot(x$lambda, xlab=xlab,ylab=ylab, ...);
    }
    if (type=="mean") {
        plot(x=x$bin.midpoints,
                       y=x$mu.x.by.bin,
                       xlab=xlab,ylab=ylab, ...);
    }
    if (type=="covariance") {
        contour(x=mids,
                        y=mids,
                        z=fit,
                        xlab="",ylab="", ... );
    }
    if (type=="correlation") {
        contour(x=mids,
                y=mids,
                z=cor(fitted(x)),
                xlab="",ylab="", ... );
    }
}