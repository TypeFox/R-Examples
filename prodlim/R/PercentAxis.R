#' Percentage-labeled axis.
#' 
#' Use percentages instead of decimals to label the an axis with a probability
#' scale .
#' 
#' 
#' @param x Side of the axis
#' @param at Positions (decimals) at which to label the axis.
#' @param \dots Given to \code{axis}.
#' @author Thomas Alexander Gerds
#' @seealso \code{\link{plot.prodlim}}
#' @keywords survival
#' @examples
#' 
#'   plot(0,0,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
#'   PercentAxis(1,at=seq(0,1,.25))
#'   PercentAxis(2,at=seq(0,1,.25))
#' 
#' @export
PercentAxis <- function(x,at,...){
  axis(x,at=at,labels=paste(100*at,"%"),...)
}
