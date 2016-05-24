#' Horizontal histogram
#' 
#' Draw a histogram with bars horizontally
#' 
#' @details Uses barplot to draw the histogram horizontally.
#
#' @return function to address y-coordinates
#' @note Doesn't work with breakpoints provided as a vector with different widths of the bars.\cr 
#'      Please do not forget to use the function for vertical positioning from the \bold{current} horizontal histogram. 
#'      If It is not working correctly, you might have the function defined from some prior horizHist result.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2011-2012
#' @seealso \code{\link{hist}}, \code{\link{barplot}}, \code{\link{axis}}
#' @keywords hplot
#' @export
#' @examples
#' 
#' # Data and basic concept
#' set.seed(8); ExampleData <- rnorm(50,8,5)+5
#' hist(ExampleData)
#' hpos <- horizHist(ExampleData) 
#' # Caution: the labels at the y-axis are not the real coordinates!
#' # abline(h=2) will draw above the second bar, not at the label value 2. 
#' # Use hpos (horizontal position), the function returned by horizHist:
#' abline(h=hpos(11), col=2, lwd=2)
#' 
#' # Further arguments
#' horizHist(ExampleData, xlim=c(-8,20)) 
#' horizHist(ExampleData, ylab="the ... argument worked!", col.axis=3) 
#' hist(ExampleData, xlim=c(-10,40)) # with xlim
#' horizHist(ExampleData, ylim=c(-10,40), border="red") # with ylim
#' hpos <- horizHist(ExampleData, breaks=20, col="orange")
#' axis(2, hpos(0:10), labels=FALSE, col=2) # another use of hpos()
#' 
#' @param Data any data that \code{\link{hist}} would take.
#' @param breaks character or numerical as explained in \code{\link{hist}}. DEFAULT: "Sturges"
#' @param freq logical. if TRUE, the histogram graphic is a representation of frequencies, the counts component of the result; 
#'        if FALSE, probability densities, component density, are plotted (so that the histogram has a total area of one). DEFAULT: TRUE
#' @param plot logical. Should histogramm be plotted? FALSE to get just the hpos function. DEFAULT: TRUE
#' @param col color. DEFAULT: par("bg")
#' @param border color of borders of bars. DEFAULT: par("fg")
#' @param las integer. Label axis style. DEFAULT: 1
#' @param xlab character. Label for x-axis. DEFAULT: "absolute frequency"
#' @param main character. Title for graphic. DEFAULT: "Histogram of substitute(Data)"
#' @param ylim numerical vector of two elements. Y-axis limits. DEFAULT: range of data
#' @param labelat numerical vector. Position of Y-Axis labels. DEFAULT: pretty(ylim)
#' @param labels numerical or character. The labels themselves. DEFAULT: labelat
#' @param \dots further arguments passed to \code{\link{barplot}} and \code{\link{axis}}
#' 
horizHist <- function(
   Data,
   breaks="Sturges",
   freq=TRUE,
   plot=TRUE,
   col=par("bg"),
   border=par("fg"),
   las=1, 
   xlab=if(freq)"Frequency" else "Density", 
   main=paste("Histogram of",deparse(substitute(Data))),
   ylim=range(HBreaks),
   labelat=pretty(ylim),
   labels=labelat,
   ... )
{
a <- hist(Data, plot=FALSE, breaks=breaks)
HBreaks <- a$breaks
HBreak1 <- a$breaks[1]
hpos <- function(Pos) (Pos-HBreak1)*(length(HBreaks)-1)/ diff(range(HBreaks))
if(plot) 
   {
   barplot(if(freq)a$counts else a$density, space=0, horiz=TRUE, ylim=hpos(ylim), col=col, border=border, 
        xlab=xlab, main=main, ...)      
   axis(2, at=hpos(labelat), labels=labels, las=las, ...) 
   }
hpos # return(invisible(hpos))
} # End of function
