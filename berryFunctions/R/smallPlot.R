#' Inset small plot within figure
#' 
#' Inset plot with margins, background and border
#' 
#' @return parameters of small plot, invisible.
#' @section Warning: setting mai etc does not work!
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2014
#' @seealso \code{\link{colPointsHist}} for an example of usage, \code{\link[TeachingDemos]{subplot}} and \code{\link[ade4]{add.scatter}} for alternative solutions to this problem that do not set margins.
#' @keywords hplot
#' @export
#' @examples
#' 
#' # Basic usage:
#' plot(1:10)
#' smallPlot(plot(5:1) )
#' smallPlot(plot(5:1), x=c(30,80), y=30:60, bg="yellow", yaxt="n")
#' # if R warns "figure margins too large", try dragging the plot viewer bigger
#' 
#' # select focus for further add-on's:
#' points(3, 2, pch="+", cex=2)
#' smallPlot( plot(5:1), bg="blue", resetfocus=FALSE )
#' points(3, 2, pch="+", cex=2)
#' 
#' # More par settings:
#' smallPlot( plot(50:1), bg=6, mai=c(0.2, 0.3, 0.1, 0.1))
#' # If you find any more that screw things up, please let me know!
#' smallPlot( plot(5:1), bg=8, ann=FALSE)
#' smallPlot(plot(10:50)) # with default bg ("transparent"), old plot is kept
#' smallPlot(plot(10:50))
#' 
#' # complex graphics in code chunks:
#' plot(1:10)
#' smallPlot( {plot(5:1, ylab="Blubber"); lines(c(2,4,3));
#'             legend("topright", "BerryRocks!", lwd=3)    }, bg="white" )
#' 
#' # in par multiple figure, things now work as well if resetfocus=TRUE:
#' op <- par("plt")
#' par(mfrow=c(2,3))
#' for(i in 1:2) plot(cumsum(rnorm(50)))
#' smallPlot( plot(50:1), bg=6)
#' plot(3:9) # opens new window
#' smallPlot( plot(50:1), bg=6, resetfocus=FALSE)
#' points(3, 2, pch="+", cex=2)
#' plot(3:9) # plot in next window, but it is still small
#' par(plt=op)
#' plot(3:9)  # margins, las and mgp are still changed
#' 
#' 
#' @param expr expression creating a plot. Can be code within {braces}.
#' @param x,y Position of small plot, relative to current figure region (0:100). max and min from vector are taken. DEFAULT: 5-70, 50-100
#' @param x1,y1,x2,y2 Positions of topleft and bottomright corner. Replaced with x,y, kept here for backcompatibility.
#' @param mar Margin vector in relative units (0:100), thus behaves differently than \code{\link{par}(mar)}. DEFAULT: c(12, 14, 3, 3)
#' @param mgp MarGinPlacement: distance of xlab/ylab, numbers and line from plot margin, as in \code{\link{par}}, but with different defaults. DEFAULT: c(1.8, 0.8, 0)
#' @param bg Background. DEFAULT: par("bg")
#' @param border Border around inset plot. DEFAULT: par("fg")
#' @param las LabelAxisStyle. DEFAULT: 1
#' @param resetfocus reset focus to original plot? Specifies where further low level plot commands are directed to. DEFAULT: TRUE
#' @param \dots further arguments passed to \code{\link{par}. new=F} removes old plot. May mess things up - please tell me for which arguments!
#'  
smallPlot <- function(
expr,
x=c(5,70),
y=c(50,100),
x1,y1,x2,y2,
mar=c(12, 14, 3, 3),
mgp=c(1.8, 0.8, 0),
bg=par("bg"),
border=par("fg"),
las=1,
resetfocus=TRUE,
...)
{                                            #     ------------
# Input check:                               #  y1 | P1       |
if(missing(x1)) x1 <- min(x, na.rm=TRUE)     #     |          |
if(missing(x2)) x2 <- max(x, na.rm=TRUE)     #  y2 |       P2 |
if(missing(y1)) y1 <- max(y, na.rm=TRUE)     #     ------------
if(missing(y2)) y2 <- min(y, na.rm=TRUE)     #       x1    x2
# catch outside plot:
if(x1<0)  {x1 <- 0;   warning("x (",x1,") set to 0.")}
if(y2<0)  {y2 <- 0;   warning("y (",y2,") set to 0.")}
if(x2>100){x2 <- 100; warning("x (",x2,") set to 100.")}
if(y1>100){y1 <- 100; warning("y (",y1,") set to 100.")}
# control for 0:1 input:
if(diff(range(x, na.rm=TRUE)) < 1  |  diff(range(y, na.rm=TRUE)) < 1  )
   stop("x or y was probably given as coodinates between 0 and 1. They must be between 0 and 100.")
# old parameters to be restored at exit:
op <- par(no.readonly=TRUE)
# inset plot: background, border
par(plt=c(x1, x2, y2, y1)/100, new=TRUE, mgp=mgp) # plt / fig
plot.new() # code line from ade4::add.scatter
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col=bg, border=border)
# inset plot: margins
par(plt=c(x1+mar[2], x2-mar[4], y2+mar[1], y1-mar[3])/100, new=TRUE, las=las, ...)
# Actual plot:
expr
# par of small plot:
sp <- par(no.readonly=TRUE)
# par reset
if(resetfocus)
  {
  if( par("mfrow")[1]==1 & par("mfrow")[2]==1  ) par(op) # ruins multiple figure plots, so:
  else par(plt=op$plt, new=op$new, mgp=op$mgp, las=op$las)
  }
return(invisible(sp))
}

