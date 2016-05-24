#' Plot method for a survival model.
#' 
#' Plot estimated baseline survival function from an object of class
#' \code{shr}. Pointwise confidence limits are available.
#' 
#' 
#' @param x a \code{shrWeib} or a \code{shrSplines} class object (output from
#' calling \code{\link{shr}} function).
#' @param type type of function to plot. The default is "shr".
#' @param add boolean.
#' @param newdata newdata.
#' @param cause cause.
#' @param col col.
#' @param lty lty.
#' @param lwd lwd.
#' @param ylim ylim.
#' @param xlim xlim.
#' @param xlab xlab.
#' @param ylab ylab.
#' @param legend legend.
#' @param confint confint.
#' @param timeOrigin timeOrigin.
#' @param axes axes.
#' @param percent percent.
#' @param \dots other graphical parameters.
#' @return Print a plot of a suvival model.
#' @author R: Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> Fortran:
#' Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' @seealso \code{\link{plot.shr}}
#' @keywords methods
#' @examples
#' 
#' # Weibull survival model
#' library(prodlim)
#' data(testdata)
#' fit.su <- shr(Hist(time=list(l,r),id)~cov,data=testdata) 
#' 
#' # pointwise confidence limits
#' plot(fit.su)
#' 
#' # no pointwise confidence limits
#' plot(fit.su,confint=FALSE)
#' 
#'
##' @S3method plot shr
plot.shr <- function(x,type="shr",add = FALSE,newdata=NULL,cause=NULL,col,lty,lwd,
	ylim,xlim,xlab="Time",ylab,legend=TRUE,confint=TRUE,timeOrigin=0,
	axes=TRUE,percent=TRUE,...){
  	
  nlines <- 1
  Y <- list(x$surv)
  plot.times <- x$time

  # {{{  getting arguments for plot, atrisk, axes, legend, confint, marktime
  if (missing(ylab)) ylab <- switch(type,"surv"=ifelse(x$reverse==TRUE,"Censoring probability","Survival probability"),"cuminc"="Cumulative incidence","hazard"="Cumulative hazard")
  if (missing(xlab)) xlab <- "Time"
  if (missing(xlim)) xlim <- c(min(plot.times), max(plot.times))
  if (missing(ylim)) ylim <- c(0, 1)
  if (missing(lwd)) lwd <- rep(3,nlines)
  if (missing(col)) col <- 1:nlines
  if (missing(lty)) lty <- rep(1, nlines)
  if (length(lwd) < nlines) lwd <- rep(lwd, nlines)
  if (length(lty) < nlines) lty <- rep(lty, nlines)
  if (length(col) < nlines) col <- rep(col, nlines)
  axis1.DefaultArgs <- list()
  axis2.DefaultArgs <- list(at=seq(0,1,.25),side=2,las=2)
  lines.DefaultArgs <- list(type="l")
  plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,xlab = xlab,ylab = ylab)
  legend.DefaultArgs <- list(legend=names(Y),lwd=lwd,col=col,lty=lty,cex=1.5,bty="n",y.intersp=1.3,x="topright")
  confint.DefaultArgs <- list(x=x,newdata=newdata,type=type,citype="shadow",times=plot.times,cause=cause,density=55,col=col[1:nlines],lwd=rep(2,nlines),lty=rep(3,nlines))
  # }}}
  smartA <- prodlim::SmartControl(call= list(...),
                         keys=c("plot","lines","legend","confint","axis1","axis2"),
                         ignore=c("x","type","cause","newdata","add","col","lty","lwd","ylim","xlim","xlab","ylab","legend","marktime","confint","automar","atrisk","timeOrigin","percent","axes","atrisk.args","confint.args","legend.args"),
                         defaults=list("plot"=plot.DefaultArgs,
                           "lines"=lines.DefaultArgs,
                           "legend"=legend.DefaultArgs,
                           "confint"=confint.DefaultArgs,
                           "axis1"=axis1.DefaultArgs,
                           "axis2"=axis2.DefaultArgs),
                         forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),
                         ignore.case=TRUE,
                         replaceDefaults=FALSE,
                         verbose=TRUE)
  # {{{  generating an empty plot
  if (!add) {
    do.call("plot",smartA$plot)
  }
  # }}}
  # {{{  adding the lines 
  lines.type <- smartA$lines$type
  nix <- lapply(1:nlines, function(s) {
    lines(x = plot.times,y = Y[[s]],type = lines.type,col = col[s],lty = lty[s],lwd = lwd[s])
  })
  # }}}
  # {{{  axes
  if (!add) {
    if (axes){
      do.call("axis",smartA$axis1)
      if (percent & is.null(smartA$axis1$labels))
        smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
      do.call("axis",smartA$axis2)
    }
  }
  # }}}

  # {{{  pointwise confidence intervals
  i <- 1
  if (confint==TRUE) {
    switch(smartA$confint$citype,
           "shadow"={
             ccrgb=as.list(col2rgb(col,alpha=T))
             names(ccrgb) <- c("red","green","blue","alpha")
             ## ccrgb$alpha=smartA$confint$density
             ccrgb$alpha=55
             cc=do.call("rgb",c(ccrgb,list(max=255)))
             xx=x$time
             nix <- sapply(1:length(xx),function(b){
               rect(xleft=xx[b],
                    xright=xx[b+1],
                    ybottom=x$lowerSurv[b],
                    ytop=x$upperSurv[b],
                    col=cc,
                    border=NA)
             })
           },{
             lines(x=x$time,x$lower,type="s",lwd=lwd[i],col=col[i],lty=lty[i],...)
             lines(x=x$time,x$upper,type="s",lwd=lwd[i],col=col[i],lty=lty[i],...)
           })
  }
  # }}}
}
