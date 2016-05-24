#' draw smoothed lines
#' 
#' draw smoothed lines with an n-level partially transparent haze
#' 
#' @return none, draws lines
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2011/2012
#' @seealso \code{\link{lines}}, \code{\link{col2rgb}}, \code{\link{rgb}}
#' @keywords aplot
#' @export
#' @examples
#' 
#' x <- 1:5  ; y <- c(0.31, 0.45, 0.84, 0.43, 0.25)
#' plot(x,y)
#' smoothLines(x,y)
#' #png("smoothLines.png")
#' par(mar=c(2,2,2,0)+.5)
#' plot(1:100, las=1, type="n", main="usage of blines(x,y, lwd, col, n, alpha ...)")
#' abline(h=0:10*10, v=0:10*10, col=6); box()
#' for(i in 0:9) { smoothLines(x=c(0,10,25,35), y=c(i*10, i*10, i*10+12, i*10+7), lwd=i)
#'                 text(25, i*10+5, paste("n=5,lwd=", i, sep="")) }
#' for(i in 0:9) { smoothLines(x=c(40,50,65,75), y=c(i*10, i*10, i*10+12, i*10+7), n=i)
#'                 text(65, i*10+5, paste("n=",i,",lwd=1", sep="")) }
#' for(i in 0:9/20) { smoothLines(x=c(80,90,105), y=c(i*200, i*200+12, i*200+12), alpha=i)
#'                    text(90, i*200+10, paste("alpha=", i, sep=""), adj=0) }
#' text(5,10, "default", adj=c(0.5,-0.2)); text(45,50, "default", adj=c(0.5,-0.2))
#' 
#' #dev.off()
#' 
#' @param x numrical. x-coordinates. x can be a matrix, then the y coordinates are taken from the second column
#' @param y numerical. y-coordinates
#' @param lwd single integer. line width
#' @param col color. DEFAULT: 1 (black)
#' @param n single integer. number of transparent lines overlayed with sinking line widths. DEFAULT: 5
#' @param alpha Transparency of color. DEFAULT: 0.1 (very transparent)
#' @param \dots further arguments as in \code{\link{lines}}
#' 
smoothLines <- function(
x,
y,
lwd=1,
col=1,
n=5,
alpha=0.1,
...)
{
# Handling for Vector and Matrix with columns x and y
if(is.vector(x)) {x <- x; y <- y} else {y <- x[,2]; x <- x[,1]}
# plot transparent lines above each other
for(i in 1:n) {lines(x,y, col=addAlpha(col, alpha=alpha), lwd=lwd+n+1-i, ...) }
# add the original line above the rest
lines(x,y, col=col, lwd=lwd, ...)   
} # end of function
