#' Plot Function For Pick-a-Point
#'
#' Plot function for objects of class \code{"pickapoint"}.
#'
#' @param x An object of class \code{"pickapoint"}.
#' @param xlab A title for the x axis (character).
#' @param ylab A title for the y axis (character).
#' @param xlim Coordinates range for x axis (numeric vector). Determined by the range of the given data by default.
#' @param ylim Coordinates range for y axis (numeric vector). Determined by the range of the given data by default.
#' @param axlwd Axis line width (numeric vector). \code{axlwd=10} by default.
#' @param cesize Size of the conditional effect marker (numeric vector). \code{cesize=1.2} by default.
#' @param cilwd Conditional interval line width (numeric vector). \code{cilwd=5} by default.
#' @param \dots Additional arguments (not supported yet).
#'
#' @return none
#'
#' @method plot pickapoint
#' @S3method plot pickapoint
#' @examples
#' \dontrun{
#' myModel <- lm('DV ~ IV + MOD', data=someData)
#' papresults <- pickapoint(myModel, dv='DV', iv='IV', mod='MOD')
#' plot(papresults)
#' }
#' @rdname plot.pickapoint
#' @export

plot.pickapoint <- function(x, xlab='', ylab='', xlim=0, ylim=0, axlwd=10, cesize=1.2, cilwd=5, ...){

  yasfunc <- function(jn,value){
    if(jn$yas == 'ratio') return(exp(value))
    else if(jn$yas == 'prob') return(exp(value) / (1 + exp(value)))
    else if(jn$yas == 'percent') return(100*(exp(value) - 1))
    else return(value)
  }

  if(is(x,'pickapoint')) {
    if(length(x$error) > 0){
      cat('Error: ',x$error)
    } else {
      #default auto scaling
      if(missing(ylim)){
        ylim=c(floor(min(x$outcome[,6])),ceiling(max(x$outcome[,7])))
      }
      if(missing(xlim)){
        xlim=c(floor(min(x$outcome[,1])),ceiling(max(x$outcome[,1])))
      }
      plot(x=x$outcome[,1], y=x$outcome[,2], ylim=ylim, xlim=xlim, ylab="" ,xlab="",pch=15,cex=cesize)

      #default axis labels
      if(ylab==''){
        ylab=paste("Conditional Effect of", x$iv, "on", x$dv)
      }
      if(xlab==''){
        xlab=x$mod
      }

      arrows(x$outcome[,1],x$outcome[,6],x$outcome[,1],x$outcome[,7],code=3,lwd=cilwd,angle=90,col='grey70')
      axis(side=1, lwd=axlwd)
      mtext(side=1, xlab, line=2,font=2)
      axis(side=2, lwd=axlwd)
      mtext(side=2, ylab, line=2,font=2)
      abline(a=x$sline, b=0)
    }
  }
}
