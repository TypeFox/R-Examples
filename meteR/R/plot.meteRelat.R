
#' @title Plot predicted METE relationships and associated observed relationship seen in data
#'
#' @description
#' \code{plot.meteRelat} plots both the theoretical prediction and data for a 
#' \code{meteRelat} object
#'
#' @details
#' \code{plot.meteRelat} automatically extracts the prediction and data (if used 
#' in \code{meteESF}) from the \code{meteDist} object. Additional plotting 
#' arguments can be passed to \code{...}.
#' 
#' @param x a \code{meteRelat} object
#' @param th.col line color of theoretical prediction 
#' @param add.legend logical; add a legend
#' @param ... arguments to be passed to \code{plot}
# @keywords manip
#' @export
#' 
#' @examples
#' data(anbo)
#' anbo.sar <- meteSAR(anbo$spp, anbo$count, anbo$row, anbo$col, Amin=1, A0=16)
#' plot(anbo.sar)
#' 
# @return list
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso meteSAR
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.

plot.meteRelat <- function(x, add.legend=TRUE, th.col='red', ...) {
    plot(x$obs, ...)
    plot(x$pred, col=th.col, add=TRUE)
    
    if(add.legend) legend('right', c('METE prediction', 'Data'), col=c(th.col, 'black'), pch=c(NA, 1), lty=c(1, NA),bty='n')
}
