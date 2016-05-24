##' Plot cumulative residuals from a 'cumres' object
##' 
##' \code{plot} displays the observed cumulative residual process with
##' realizations under the null. 95\% prediction bands
##' 
##' 
##' @param x Object produced by the function \code{cumres}.
##' @param idx vector of numbers (or variable names) indicating which processes
##' from the \code{x} to plot.
##' @param col Color of the sample processes. By setting this parameter to "none"
##' or \code{NULL} no realizations will be drawn. The number of realizations is
##' determined by the \code{cumres}-object.
##' @param ci Type of prediction bands to plot. Defaults to none. Set to
##' \code{TRUE} to obtain simultaneous prediction bands under the null (pointwise
##' can be obtained by setting to "pointwise").
##' @param col.ci Color of prediction band.
##' @param col.alpha Degree of transparency (0-1) of the prediction bands.
##' @param lty.ci Line type of prediction band.
##' @param level The required prediction level.
##' @param legend Type of legend where "type1" gives p-values of GoF-tests and
##' "type2" gives usual type of legends.
##' @param xlab Optional label of x-axis
##' @param ylab Optional label of y-axis
##' @param vs Label of predictor
##' @param ylim Range of y axis
##' @param title Main title
##' @param ... Additional arguments passed to the plot-routine.
##' @author Klaus K. Holst <kkho@@biostat.ku.dk>
##' @seealso \code{\link[gof]{cumres}}
##' @keywords hplot regression
##' @examples
##' 
##' n <- 500; x <- abs(rnorm(n,sd=0.2))+0.01; y <- sqrt(x) + rnorm(n,sd=0.2)
##' l <- lm(y ~ x)
##' g <- cumres(l, R=500)
##' plot(g, idx=1, ci="sim", col=NULL, col.ci="purple", legend="type2")
##' @method plot cumres
##' @export
plot.cumres <- function(x, idx=1:length(x$variable),
                        col=c("grey"),
                        ci=TRUE,
                        col.ci="darkblue", col.alpha=0.3, lty.ci=0, level=0.95,
                        legend=c("type1","type2","none"), xlab, ylab,
                        vs=TRUE,
                        ylim=NULL,
                        title,
                        ...) {
  ylim. <- ylim
  newylab <- missing(ylab)
  newxlab <- missing(xlab)
  for (i in idx) {
    legendtxt <- c(); legendpch <- c(); legendcol <- c(); legendlty <- c(); legendlwd <- c(); legendcex <- c()
    if (is.null(ylim)) {
      ylim. <- max(abs(range(x$W[,i])))*2*c(-1,1)
    }
    ## Observed process
    main <- ""
    if (newxlab) {
      xlab <- x$variable[i]; 
      if (x$type[i]=="score") {
        main <- xlab; xlab <- "Time";
      }
    }
    if (newylab) {
      ylab <- substitute(expression(W[p](x)),list(p=x$variable[i]))
    }
    legendtxt <- c(legendtxt, "Observed"); legendpch <- c(legendpch,-1); legendcol <- c(legendcol,1); legendlty <- c(legendlty,1); legendlwd <- c(legendlwd,2); legendcex <- c(legendcex,1);
    x0 <- x$x[,i]
    if (!vs) {
      x0 <- 1:length(x0)
      xlab <- "Observation"
    }

    
    sampleproc <- function() {
      ## Sample processes
      if (col!="none" && !is.null(col)) {
        legendtxt <- c(legendtxt, "MC sample"); legendpch <- c(legendpch,-1); legendcol <- c(legendcol,col); legendlty <- c(legendlty,1); legendlwd <- c(legendlwd,1); legendcex <- c(legendcex,1);
        for (k in 1:ncol(x$What[[i]])) {
          lines(x$What[[i]][,k] ~ x0, type="s", col=col, lwd=1)
        }; lines(x$W[,i] ~ x0, type="s", lwd=2)
      }
    }
    predband <- function() {
      ## Prediction bandds
      if ( (ci[1]!="none" && !is.null(ci) && ci[1]!=0) || (ci==TRUE) ) {
        ##    if ((ci[1]!="none" && !is.null(ci))) {
        if (col.alpha==0)
          col.trans <- col.ci
        else 
          col.trans <- sapply(col.ci, FUN=function(x) do.call(rgb,as.list(c(col2rgb(x)/255,col.alpha))))
        if (ci[1]=="pointwise")
          myCI <- confint(x,parm=i,cval=qnorm(1-(1-level)/2))
        else
          myCI <- confint(x,parm=i,level=level)      
        mystepf <- with(myCI, stepfun(t,c(0,yu)));
        t <- c();
        epsilon <- 1e-9
        for (k in 1:length(myCI$t)) {
          t <- c(t, myCI$t[k]-epsilon, myCI$t[k])
        }; t <- t[-1]
        yu <- mystepf(t)
        legendtxt <- c(legendtxt, "95% prediction band"); legendpch <- c(legendpch,15); legendcol <- c(legendcol,col.trans); legendlty <- c(legendlty,0); legendlwd <- c(legendlwd,0); legendcex <- c(legendcex,2);
        
        lines(yu ~ t, lwd=1, col=col.ci, lty=lty.ci)
        lines(-yu ~ t, lwd=1, col=col.ci, lty=lty.ci)
        tt <- c(t, rev(t))
        yy <- c(yu, rev(-yu))
        polygon(tt,yy, col=col.trans, lty=0)      
        ##as.list(c(col2rgb("darkblue"),10)/255))),
      }
    }

    with(x, plot(W[,i] ~ x0, type="n", lwd=2, ylab=ylab, ylim=ylim.,xlab=xlab,main=main));    
    if (col.alpha==0) {
      with(x, lines(W[,i] ~ x0, type="s", lwd=2));
      predband()
      sampleproc()
    } else {
      with(x, lines(W[,i] ~ x0, type="s", lwd=2));
      sampleproc()      
      predband()
    }
    
    if (!missing(title)) {
      graphics::title(title)
    } else {
      if (!is.null(x$response))
        graphics::title(x$response)
    }
    
    if (!is.null(legend) && legend[1]!="none" && (legend!=F)) {
      if (legend[1]=="type1")
        legend("topright", c(paste("KS-test: p=",x$KS[i],sep=""),paste("CvM-test: p=",x$CvM[i],sep="")), bg="white")
      else
        legend("topright", legendtxt, lty=legendlty, pch=legendpch, col=legendcol, lwd=legendlwd, pt.cex=legendcex, bg="white")
    }
    ylim. <- NULL
  }
  invisible(x)
}



