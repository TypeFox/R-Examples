# pagenum.r

#' Put Page Numbers on Graphics
#'
#' @details
#' Use \code{setPagenum} to create a global variable with the page number.
#'
#' If \code{pagenum()} is called without an argument, the value of
#' \code{options()$pagenum} is used to determine the page number.
#'
#' Each time \code{pagenum()} is called, \code{options()$pagenum} is
#' automatically incremented by 1.
#'
#' @rdname pagenum
#' @aliases pagenum setPagenum getPagenum
#'
#' @param num The number to put on the page.  If no number is given,
#' the value of \code{options()$pagenum} is used.
#' @param text The text to use in front of the page number.
#' @param date If FALSE (default), do not add a date below the page number.
#' @param date.format The format to use for the date.
#' @param x Horizontal position of timestamp, in [0,1]. Default .03
#' @param y Vertical position of timestamp, in [0,1]. Default .03
#' @param just Jufstification.  Default c('left','bottom')
#' @param col Color to use for the text.
#' @param cex Character expansion. Default 0.75.
#'
#' @return Returns the value stored by the global variable.
#' @author Kevin Wright
#' @references
#' Mark Heckmann (2009).
#' R: Good practice - adding footnotes to graphics.
#' \url{https://ryouready.wordpress.com/2009/02/17/r-good-practice-adding-footnotes-to-graphics/}
#'
#' @examples
#' # base graphics
#' setPagenum(1)
#' plot(Sepal.Length~Sepal.Width, data=iris, col=Species, pch=19)
#' pagenum()
#'
#' # lattice, date
#' setPagenum(getPagenum()+1) # Manual increment
#' require(lattice)
#' xyplot(Sepal.Length~Sepal.Width, data=iris, groups=Species)
#' pagenum(date=TRUE)
#'
#' # ggplot2, top-right
#' require(ggplot2)
#' ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length,
#'   color=Species)) + geom_point() + theme_classic()
#' pagenum(text="ABC Corp - ", date=TRUE,
#' x=.95, y=.95, just=c('right','top'))
#'
#' # multiple figures, 'draft' watermark
#' op = par(mfrow=c(1,2))
#' plot(Sepal.Length ~ Sepal.Width, data=iris, col=Species, pch=19)
#' plot(Petal.Length ~ Petal.Width, data=iris, col=Species, pch=19)
#' par(op)
#' pagenum(num="", text="Figures 2a, 2b")
#' pagenum(num="", text="Draft",
#'         x=.5, y=.95, just=c('center','top'),
#'         col="wheat", cex=3)
#'
#' @import grid
#' @export
pagenum <- function(num, text="Page", date=FALSE,
                    date.format,
                    x=.03, y=.03, 
                    just=c("left","bottom"),
                    col="gray50", cex=.75){

  # Can't get roxygen and Rd_parse to both work on "%d %b %Y %H:%M:%S"
  # so set it manually here instead of as a default in the function call.
  if(missing(date.format))
    date.format <- "%d %b %Y %H:%M:%S"

  # If num is missing, get value from options()
  if(missing(num)) num <- getPagenum()

  pn <- paste(text,num)

  if(date)
    pn <- paste(pn, "\n", format(Sys.time(), date.format), sep="")

  gp=gpar(cex=cex, col=col)

  # Need to clip if there are multiple base figures
  # See: http://tolstoy.newcastle.edu.au/R/help/06/06/30031.html
  grid.clip() 
  pushViewport(viewport(x, y, width=stringWidth(pn),
                        height=unit(2,"lines"),
                        name="pagenum", gp=gp))
  grid.text(pn, gp=gp, just=just)
  popViewport()
  
  # If num is numeric, increment page number counter for next time.
  if(is.numeric(num)) setPagenum(num+1)

  invisible()
}

#' @rdname pagenum
#' @export
setPagenum <- function(num=1){
  options(pagenum=num)
  invisible()
}

#' @rdname pagenum
#' @return Returns the value of options()$pagenum
#' @export
getPagenum <- function(){
  # If the setting is missing, create it and set it to 1
  if(is.null(options()$pagenum)) options(pagenum=1)

  return(options()$pagenum)
}
# ----------------------------------------------------------------------------

if(FALSE){

  # Here is an example of annotating EACH figure of multi-figure lattice graphics
  # http://r.789695.n4.nabble.com/adding-text-to-the-top-corner-of-a-lattice-plot-td3092079.html
  require("grid")
  require("lattice")
  require("gridExtra") # for grid.arrange
  pgfun <- function(mark) {
    function(n) grid.text(label = mark,
                          x = unit(0.1, "npc"), y = unit(0.9, "npc"))
  }
  # The 'page' argument must be a function with 1 argument, page number.
  # Some trickery with environments to pass in the text string.
  grid.arrange(nrow = 2,
               xyplot(demand ~ Time, BOD, page = pgfun("A")),
               xyplot(demand ~ Time, BOD, page = pgfun("B")),
               xyplot(demand ~ Time, BOD, page = pgfun("C"))
             )

}
