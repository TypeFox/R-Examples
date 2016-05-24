#' zoom in originally static x11 graphics
#' 
#' zoom in x11 graphics - uses locator to define region to zoom into
#' 
#' @return none, works in existing graphics
#' @note Extensive testing is yet to be done!
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, ca 2012
#' @seealso shapeZoom in \url{https://github.com/brry/shapeInteractive}, \code{\link{colPoints}}, \code{\link{locator}}
#' @keywords aplot iplot dynamic
#' @export
#' @examples
#' 
#' ## Examples rely on locator, so can't be checked in non-interactive R use.
#' \dontrun{
#' ## Rcmd check --as-cran doesn't like to open external devices,
#' ## so this example is excluded from running in the checks.
#' a <- rnorm(90); b <- rexp(90)
#' windows(record=TRUE) # turn recording on
#' plot(a,b, las=1)
#' pointZoom(a,b, col=2, expr="abline(v=0)")
#' # now scroll through the plots (Pg Up and Pg Dn) to unzoom again.
#' 
#' d <- data.frame(a,b)
#' class(d)
#' plot(d)
#' pointZoom(d)
#' }
#' 
#' @param x same x coordinates as in current plot. x can be a matrix, then the y (and z) coordinates are taken from the second (and third) column.
#' @param y ditto
#' @param z if using colpoints, z-value
#' @param Time Duration of zooming (speed) in seconds. DEFAULT: 1
#' @param steps number of single zoomlevels. DEFAULT: 30
#' @param las label axis style, see \code{\link{par}}. DEFAULT: 1
#' @param usecolp logical: use \code{\link{colPoints}} when zooming? DEFAULT: FALSE
#' @param xlab xlabel See \code{\link{plot}}. DEFAULT: substitute(x)
#' @param ylab dito
#' @param quiet logical. Should notifications (instructions) be written to the console? DEFAULT: FALSE
#' @param expr Characterized Expression to be executed after each plot, eg. \code{expr='abline(h=3)'}
#' @param \dots further arguments passed to \code{\link{plot}} or \code{\link{colPoints}}.
#' 
pointZoom <- function(
x,
y=NA,
z=NA,
Time=1,
steps=30,
las=1,
usecolp=FALSE,
xlab=substitute(x),
ylab=substitute(y),
quiet=FALSE,
expr,
...)
{
if(interactive()){ # to silence the R CMD check warnings
if(!quiet){
  legend("top", "Instructions appear in console", bty="n", text.col="orange")
  message("Please select the area to zoom to in the graphics window.")
  message("first klick topleft, then bottomright."); flush.console()   
  } # if notify end
w <- locator(2)
u <- par()$usr
if(w$x[1] > w$x[2] | w$y[1] < w$y[2])
  {
  message("wrong selection!")
  message("first klick topleft, then bottomright of the area to zoom to.")
  flush.console(); w <- locator(2)
  }
# if x is matrix:
if(inherits(x,c("matrix", "data.frame", "array")) )
   {y <- x[,2]; if(ncol(x)>2) z <- x[,3]; x <- x[,1]}
#
if (usecolp)
  {polygon(c(w$x, rev(w$x)), rep(w$y, each = 2))
  Sys.sleep(1)
  colPoints(x, y, z, add=FALSE, xlim=w$x, ylim=rev(w$y), las=las, ylab=ylab, xlab=xlab, ...)
  } else  {
X1 <- c(u[1]+(w$x[1]-u[1])*1:steps/steps)
X2 <- c(u[2]-(u[2]-w$x[2])*1:steps/steps)
Y1 <- c(u[3]+(w$y[2]-u[3])*1:steps/steps)
Y2 <- c(u[4]-(u[4]-w$y[1])*1:steps/steps)
for ( i in 1:steps) 
   {
   rect(xleft=w$x[1], ybottom=w$y[1], xright=w$x[2], ytop=w$y[2])
   Sys.sleep(Time/steps)
   plot(x, y, xlim=c(X1[i], X2[i]), ylim=c(Y1[i], Y2[i]), las=las,
   ylab=ylab, xlab=xlab ,  yaxs="i", xaxs="i", ...)
   rect(xleft=w$x[1], ybottom=w$y[1], xright=w$x[2], ytop=w$y[2])
   if(!missing(expr)) eval(parse(text=expr))
   } # loop end
}
# if(!quiet) message("Tell me if this was helpful: berry-b@gmx.de") 
} # end if interactive
} # function end
