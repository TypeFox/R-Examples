#' Sequential color palette
#' 
#' Sequential color palette from yelow to red or yellow to blue or custom colors.
#' 
#' @return Character string vector with color names
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Jan 2016
#' @seealso \code{\link{showPal}}, \code{\link{divPal}}, \code{\link{addAlpha}}, 
#'          \code{\link{colorRampPalette}}, package \code{RColorBrewer}
#' @keywords color dplot
#' @export
#' @examples
#' 
#' plot(rep(1,12),  pch=16, cex=5, col=seqPal(12), xaxt="n")
#' showPal()
#'  
#' @param n Number of colors. DEFAULT: 12
#' @param reverse Reverse colors? DEFAULT: FALSE
#' @param alpha Transparency (0=transparent, 1=fully colored). DEFAULT: 1
#' @param extr Should colors span possible range more extremely? 
#'       If TRUE, it has very light yellow and very dark blue values included, 
#'       using the result from \code{RColorBrewer::brewer.pal(9, "YlGnBu")}. DEFAULT: FALSE
#' @param yb Should colors be in yellow-blue instead of the internal (nice) default? DEFAULT: FALSE
#' @param yr Should colors be in yellow-red instead of the default? DEFAULT: FALSE
#' @param colors If not NULL, a color vector used in \code{\link{colorRampPalette}}. DEFAULT: NULL
#' @param \dots Further arguments passed to \code{\link{colorRamp}}
#' 
seqPal <- function(
n=12,
reverse=FALSE,
alpha=1,
extr=FALSE,
yb=FALSE,
yr=FALSE,
colors=NULL,
...
)
{
cols <- c("#FFFFC6","#CAE9AE","#85CDBA","#4DB6C6","#327EBD","#22329A")
if(extr) cols <- c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4",
                    "#1D91C0", "#225EA8", "#253494", "#081D58")
if(yb) cols <- c("yellow","blue")
if(yr) cols <- c("yellow","red")
if(!is.null(colors)) cols <- colors
if(reverse) cols <- rev(cols)
outcols <- colorRampPalette(cols)(n)
if(alpha!=1) outcols <- addAlpha(outcols, alpha)
outcols
}
