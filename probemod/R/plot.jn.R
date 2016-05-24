#' Plot Function For Johnson-Neyman
#'
#' Plot function for objects of class \code{"jn"}.
#'
#' @param x An object of class \code{"jn"}.
#' @param xlab A title for the x axis (character).
#' @param ylab A title for the y axis (character).
#' @param xlim Coordinate range for x axis (numeric vector). Determined by the range of the given data by default.
#' @param ylim Coordinate range for y axis (numeric vector). Determined by the range of the given data by default.
#' @param axlwd Axis line width (numeric vector). \code{axlwd=10} by default.
#' @param celwd Conditional effect line width (numeric vector). \code{celwd=10} by default.
#' @param cblwd Conditional band line width (numeric vector). \code{cblwd=8} by default.
#' @param \dots Additional arguments (not supported yet).
#'
#' @return none
#'
#' @method plot jn
#' @S3method plot jn
#' @examples
#' \dontrun{
#' myModel <- lm('DV ~ IV + MOD', data=someData)
#' jnresults <- jn(myModel, dv='DV', iv='IV', mod='MOD')
#' plot(jnresults)
#' }
#' @rdname plot.jn
#' @export

plot.jn <- function(x, xlab='', ylab='', xlim=0, ylim=0, axlwd=10, celwd=10, cblwd=8, ...){

  if(is(x,'jn')) {
    if(length(x$error) > 0){
      cat('Error: ',x$error)
    } else {
      #default auto scaling
      if(missing(ylim)){
        ylim=c(floor(min(x$plot$data$llci)),ceiling(max(x$plot$data$ulci)))
      }
      if(missing(xlim)){
        xlim=c(floor(min(x$plot$data$x)),ceiling(max(x$plot$data$x)))
      }
      plot(x=x$plot$data$x, y=x$plot$data$y, ylim=ylim, xlim=xlim, ylab="" ,xlab="")

      #default axis labels
      if(ylab==''){
        ylab=paste("Conditional Effect of", x$iv, "on", x$dv)
      }
      if(xlab==''){
        xlab=x$mod
      }

      cat(paste('Values of', x$mod, 'indicated by the shaded region\n'))
      srvect <- vector()
      srvect = rbind(srvect, unlist(lapply(x$plot$summary$lower, '[[', 1)))
      srvect = rbind(srvect, unlist(lapply(x$plot$summary$upper, '[[', 1)))
      rownames(srvect) <- c('Lower Bound:','Upper Bound:')
      #srvect[,'y'] <- yasfunc(x,srvect[,'y'])
      print(srvect)
      for(i in 1:nrow(x$plot$signintervals)) {
        ypolycoords <- x$plot$data$x[x$plot$signintervals[i,1]:x$plot$signintervals[i,2]]
        polygon(c(ypolycoords,rev(ypolycoords)),c(x$plot$data$llci[x$plot$signintervals[i,1]:x$plot$signintervals[i,2]],x$plot$data$ulci[x$plot$signintervals[i,2]:x$plot$signintervals[i,1]]),col = "grey50", density = c(10, 20), angle = c(-45, 45), border = FALSE)
      }
      lines(x=x$plot$data$x,y=x$plot$data$y,lwd=celwd)
      lines(x=x$plot$data$x, y=x$plot$data$ulci, col="grey50",lty=2,lwd=cblwd)
      lines(x=x$plot$data$x, y=x$plot$data$llci, col="grey50",lty=2,lwd=cblwd)
      axis(side=1, lwd=axlwd)
      mtext(side=1, xlab, line=2,font=2)
      axis(side=2, lwd=axlwd)
      mtext(side=2, ylab, line=2,font=2)
      abline(a=x$plot$sline, b=0)
    }
  }
}
