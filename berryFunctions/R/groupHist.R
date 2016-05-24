#' Histogramm for classes
#' 
#' Improvement of \code{tapply(x, g, hist)} with x and g taken from a data.frame
#' 
#' @details Uses \code{\link{split}} to categorize into groups.
#' 
#' @return NULL, used for plotting
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Jan 2015
#' @seealso \code{\link{hist}}, \code{\link{tapply}}
#' @keywords dplot hplot distribution
#' @export
#' @examples
#' 
#' groupHist(chickwts, weight, "feed", col=2)
#' groupHist(chickwts, "weight", "feed", col=2, unit="grams at age 6 weeks")
#' groupHist(chickwts, weight, feed, col=2, breaks=20, main="Hi there")
#' groupHist(iris, Petal.Width, Species)
#' 
#' @param d data.frame object name
#' @param x column name of variable of interest
#' @param g column name of groups (\code{INDEX in \link{tapply}, f in \link{split}})
#' @param xlab,ylab axis labels. DEFAULT: ""
#' @param las LabelAxisStyle, see \code{\link{par}}. DEFAULT: 1, means numbers on y-axis upright
#' @param main Main title, internal default based on \code{d, x, unit} and \code{g}. DEFAULT: NULL
#' @param unit Unit to be written into the default title. DEFAULT: NA
#' @param \dots further arguments passed to \code{\link{hist}}
#' 
groupHist <- function(
d,
x,
g,
xlab="", ylab="",
las=1,
main=NULL,
unit=NA,
...
)
{
# Input can be quoted, must not:
x <- gsub("\"", "", deparse(substitute(x)), fixed=TRUE)
g <- gsub("\"", "", deparse(substitute(g)), fixed=TRUE)
# get break values from full dataset,
# ignore messages like: argument 'col' is not made use of
suppressWarnings( h <- hist(d[,x], plot=FALSE, ...)    )
# prepare dimensions:
names <- unique(d[,g])
n <- length(names)
op <- par(mfrow=panelDim(n, landscape=TRUE), mar=c(3,3,2,1), mgp=c(3,0.7,0),
          oma=c(0,0,2,0))
on.exit(par(op))
# perform grouping:
dd <- split(d[,x], d[,g])
# main title:
unit2 <- if(is.na(unit)) "" else paste0("  [", unit, "] ")
if(is.null(main))
   main <- paste0("Histograms of ",x, unit2," in ",substitute(d),", grouped by ", g)
# actual plotting:
for(i in 1:n)
  {
  defaults <- list(x=dd[[i]], main=names(dd)[i], breaks=h$breaks, xlab=xlab, ylab=ylab, las=las)
  do.call(hist, owa(defaults, list(...), "breaks"))
#  hist(x=dd[[i]], main=names(dd)[i], breaks=h$breaks, xlab=xlab, ylab=ylab, las=las, ...)
  box("figure")
  }
title(main=main, outer=TRUE, line=0.7)
}
