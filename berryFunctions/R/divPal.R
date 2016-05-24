#' Diverging color palette
#' 
#' Diverging color palette: brown to blue, light colors in the middle, darker
#' at the extremes, good for displaying values in two directions
#' 
#' @return Character string vector with color names
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Jan 2016
#' @seealso \code{\link{showPal}}, \code{\link{seqPal}}, \code{\link{addAlpha}}, 
#'          \code{\link{colorRampPalette}}, package \code{RColorBrewer}
#' @references Originally in 12 shades in the IPCC Assesment Report 5 Chapter 12 Fig 12.22, \url{http://www.ipcc.ch/report/ar5/wg1/}
#' @keywords color dplot
#' @export
#' @examples
#'
#' plot(rep(1,12), pch=16, cex=5, col=divPal(12), xaxt="n")
#' showPal()
#' 
#' @param n Number of colors. DEFAULT: 12
#' @param reverse Reverse colors? DEFAULT: FALSE
#' @param alpha Transparency (0=transparent, 1=fully colored). DEFAULT: 1
#' @param rwb Should colors be in red-white-blue instead of brown-blue? DEFAULT: FALSE
#' @param ryb Use red-yellow-blue instead of the default, with "khaki" in the center. DEFAULT: FALSE
#' @param colors If not NULL, a color vector used in \code{\link{colorRampPalette}}. DEFAULT: NULL
#' @param \dots Further arguments passed to \code{\link{colorRamp}}
#' 
divPal <- function(
n=12,
reverse=FALSE,
alpha=1,
rwb=FALSE,
ryb=FALSE,
colors=NULL,
...
)
{
cols <- c("#9B5523", "#B16A32", "#CA8448", "#F4C882", "#F1DB99", "#FBF5B4",
          "#C3E2C0", "#96D1A7", "#46BEA0", "#4984A0", "#4984A0", "#0B3A5B")
if(rwb) cols <- c("red","white","blue")
if(ryb) cols <- c("red","khaki1","blue")
if(!is.null(colors)) cols <- colors
if(reverse) cols <- rev(cols)
outcols <- colorRampPalette(cols)(n)
if(alpha!=1) outcols <- addAlpha(outcols, alpha)
outcols
}
