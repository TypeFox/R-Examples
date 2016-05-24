#' Legend for colPoints
#' 
#' Adds legends to plots created or enhanced with \code{\link{colPoints}}
#' 
#' @note \code{x1,x2,y1,y2,labelpos,titlepos,title} have different defaults when \code{horizontal=FALSE}
#' @return invisible list of par of \code{\link{smallPlot}}, adds legend bar to current plot
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2012-2014
#' @seealso \code{\link{colPoints}} (section examples) for real life example
#' @keywords aplot color
#' @export
#' @examples
#' 
#' z <- rnorm(50)
#' plot(1:10)
#' colPointsLegend(z=z)
#' colPointsLegend(z=z, titlepos=2)
#' colPointsLegend(z=z, horiz=FALSE) # note the different defaults
#' # positioning relative to plot:
#' colPointsLegend(z=z, x1=5, x2=30, y1=90, y2=70, title="Booh!", density=FALSE)
#' # Denote values outside of Range wit a triangle:
#' colPointsLegend(z=z, Range=c(-1,3), x1=20, y1=60, y2=40, triangle=c(0,0.5))
#' colPointsLegend(z=z, horiz=FALSE, x1=70, y1=60, plottriangle=TRUE, density=FALSE)
#' ?colPoints # example section for actual usage
#' 
#' @param z Values of third dimension used in \code{\link{colPoints}}, can be matrix, vector, etc, but must be numeric
#' @param Range Ends of color bar for method=equalinterval. DEFAULT: range(z, finite=TRUE)
#' @param nbins Number of classes (thus, colors). DEFAULT: 40
#' @param colors Color vector. DEFAULT: \code{\link{rainbow}} from blue (lowest) to red (highest value in Range)
#' @param bb Borders of bins for the legend (key). DEFAULT: seqR(Range, length.out=nbins+1)
#' @param at Positions of legend labels. DEFAULT: pretty2(Range)
#' @param labels Labels that are written at the positions of \code{at}. DEFAULT: at
#' @param adj label adjustment parallel to legend bar (only one number!). DEFAULT: 0.5
#' @param bg Background behind key, labels and title. DEFAULT: "white"
#' @param x1,y1 Topleft relative coordinates (0:100) of inset plot, see \code{\link{smallPlot}}. DEFAULT: 60,99
#' @param x2,y2 Bottomright -"-. DEFAULT: 98,88
#' @param mar Margins for \code{\link{smallPlot}} in relative values (0:100). DEFAULT: internal calculations based on title, labelpos and titlepos.
#' @param mgp MarGinPlacement: distance of xlab/ylab, numbers and line from plot margin, as in \code{\link{par}}, but with different defaults. DEFAULT: c(1.8, 0.6, 0)
#' @param sborder Border around inset subplot. DEFAULT: NA
#' @param resetfocus Reset focus to original plot? Specifies where further low level plot commands are directed to. DEFAULT: TRUE
#' @param plottriangle Should triangles be plotted at the end of the legend for values outside Range? TRUE if missing but triangle is given. DEFAULT: FALSE
#' @param triangle Percentage of bar length at lower and upper end for triangles (can be a vector with two different values). DEFAULT: 0.14
#' @param tricol Triangle colors for lower and upper end. DEFAULT: c(1,8)
#' @param density Plot kernel density line? arguments passed to \code{\link{density}}. DEFAULT: NULL
#' @param lines Plot black lines in the color bar at \code{at}? DEFAULT: TRUE
#' @param atminmax Should the extrema of the legend be added to \code{at}? DEFAULT: FALSE
#' @param horizontal Horizontal bar? if FALSE, a vertical bar is drawn. DEFAULT: TRUE
#' @param labelpos Position of labels relative to the bar. Possible: 1 (below), 2 (left), 3 (above), 4 (right), 5(on top of bar). DEFAULT: 1
#' @param titlepos Position of title -"-. DEFAULT: 3
#' @param title Legend title. DEFAULT: "Legend"
#' @param las LabelAxisStyle. DEFAULT: 1
#' @param \dots Further arguments passed to \code{\link{text}} and \code{\link{strwidth}}, e.g. cex, srt, font, col. But NOT adj!
#' 
colPointsLegend <- function(
z,
Range=range(z, finite=TRUE),
nbins=40,
colors=seqPal(nbins),
bb=seqR(Range, length.out=nbins+1),
at=pretty2(Range),
labels=at,
adj=0.5,

bg="white",
x1=60,
y1=99,
x2=x1+38,
y2=y1-11,
mar,
mgp=c(1.8, 0.6, 0),
sborder=NA,
resetfocus=TRUE,

plottriangle=FALSE,
triangle=0.14,
tricol=c(1,8),
density=NULL,
lines=TRUE,
atminmax=FALSE,
horizontal=TRUE,
labelpos=1,
titlepos=3,
title="Legend",
las=1,
...)
{
# ------------------------------------------------------------------------------
z <- as.numeric(z)
# input checks:
if(any(diff(bb)<0)) stop("Breaks 'bb' (bin borders) have to be in ascending order.")
if(missing(nbins) & !missing(colors)) nbins <- length(colors)
if(length(colors) != nbins) stop("Number of colors is not equal to number of classes.")
# extend labels and at:
if(atminmax) labels <- c( signif(head(bb,1),2), labels, signif(tail(bb,1),2) ) ### & length(labels)!=length(at)
if(atminmax) at <- c( head(bb,1), at, tail(bb,1) )
if(length(labels)!=length(at)) stop("labels and at do not have the same length")
# vertical default placement:
if(!horizontal){
if(missing(x1)) x1 <- 88
if(missing(y1)) y1 <- 70
if(missing(x2)) x2 <- pmin(x1+11, 100)
if(missing(y2)) y2 <- pmax(y1-40, 0)
if(missing(labelpos)) labelpos <- 2
if(missing(titlepos)) titlepos <- 3
if(missing(title)) title <- "Key"
}
# triangle preparation:
if(!missing(triangle) & missing(plottriangle)) plottriangle <- TRUE
if(plottriangle)
  {
  if(!is.numeric(triangle)) stop("triangle must be numeric.")
  if(any(triangle>2 | triangle<0)) stop("Values in triangle must be between 0 and 2")
  triangle <- rep(triangle, length.out=2)
  tricol   <- rep(tricol  , length.out=2)
  # coordinates of triangle points
  barlength <- (tail(bb,1) - bb[1])
  trimin <- bb[1] - barlength*triangle[1]
  trimax <- tail(bb,1) + barlength*triangle[2]
  plotrange <- c(trimin, trimax)
  }  else
  plotrange <- c(bb[1], tail(bb,1))
# margin preparation:
if(missing(mar))
{
mar <- c(0,0,0,0)
wt <- 1.4*100*max( strwidth(c(labels, title), units="figure", ...))
ht <- 1.5*100*max(strheight(c(labels, title), units="figure", ...))
if(labelpos==2 | titlepos==2) mar[2] <- wt
if(labelpos==4 | titlepos==4) mar[4] <- wt
if(labelpos==1 | titlepos==1) mar[1] <- ht
if(labelpos==3 | titlepos==3) mar[3] <- ht
} # if mar is specified, it is used, of course.
# subplot setup:
smallPlot(x1=x1, y1=y1, x2=x2, y2=y2, mar=mar, mgp=mgp, bg=bg,
  border=sborder, las=las, resetfocus=resetfocus, expr={
if(horizontal) # ---------------------------------------------------------------
  {
  plot.window(xlim=plotrange, ylim=c(0,1), xaxs="i", yaxs="i")
  # actually plot legend color fields:
  for(i in 1:length(colors))
    rect(xleft=bb[i], xright=bb[i+1], ybottom=0, ytop=1, col=colors[i], border=NA)
  # triangle:
  if(plottriangle)
     {
     polygon(c(bb[1],bb[1],trimin),       c(0,1,0.5), col=tricol[1], border=NA)
     polygon(c(rep(tail(bb,1),2),trimax), c(0,1,0.5), col=tricol[2], border=NA)
     }
  # lines
  if(lines) segments(x0=at, y0=0, y1=1)
  # prepare label adjustment:
  if(labelpos==1) { y <- -0.1 ; vadj <- 1   } else
  if(labelpos==3) { y <-  1.1 ; vadj <- 0   } else
  if(labelpos==5) { y <-  0.5 ; vadj <- 0.5 } else
  stop("Wrong labelpos. Possible in horizontal legend: 1 (below legend bar), 3 (above), and 5 (on top).")
  # actually write labels:
  text(x=at, y=y, labels=labels, adj=c(adj, vadj), xpd=TRUE, ...)
  # prepare title adjustment:
  pu <- par("usr")[1:2]
  if(titlepos==1) {x <- mean(pu); y <- -0.2; hadj <- 0.5; vadj <- 1   } else
  if(titlepos==2) {x <-    pu[1]; y <-  0.5; hadj <- 1  ; vadj <- 0.5 } else
  if(titlepos==3) {x <- mean(pu); y <-  1.2; hadj <- 0.5; vadj <- 0   } else
  if(titlepos==4) {x <-    pu[2]; y <-  0.5; hadj <- 0  ; vadj <- 0.5 } else
  if(titlepos==5) {x <- mean(pu); y <-  0.5; hadj <- 0.5; vadj <- 0.5 } else
  stop("Wrong titlepos. Must be integer between 1 and 5.")
  # actually write title:
  text(x=x, y=y, labels=title, adj=c(hadj, vadj), xpd=TRUE, ...)
  # kernel density:
  if(is.list(density) | is.null(density) | isTRUE(density) )
    {
    dp <- do.call(stats::density, args=owa(list(x=z, na.rm=TRUE), density))
    lines(dp$x, dp$y/max(dp$y))
    }
  }
else # if not horizontal, thus if vertical -------------------------------------
  {
  plot.window(ylim=plotrange, xlim=c(0,1), xaxs="i", yaxs="i")
  # actually plot legend color fields:
  for(i in 1:length(colors))
    rect(ybottom=bb[i], ytop=bb[i+1], xleft=0, xright=1, col=colors[i], border=NA)
  # triangle:
  if(plottriangle)
     {
     polygon(c(0,1,0.5), c(bb[1],bb[1],trimin),       col=tricol[1], border=NA)
     polygon(c(0,1,0.5), c(rep(tail(bb,1),2),trimax), col=tricol[2], border=NA)
     }
  # lines
  if(lines) segments(y0=at, x0=0, x1=1)
  # prepare label adjustment:
  if(labelpos==2) { x <- -0.1 ; hadj <- 1   } else
  if(labelpos==4) { x <-  1.1 ; hadj <- 0   } else
  if(labelpos==5) { x <-  0.5 ; hadj <- 0.5 } else
  stop("Wrong labelpos. Possible in vertical legend: 2 (left of legend bar), 4 (right), and 5 (on top).")
  # actually write labels:
  text(x=x, y=at, labels=labels, adj=c(hadj, adj), xpd=TRUE, ...)
  # prepare title adjustment:
  pu <- par("usr")[3:4]
  if(titlepos==1) {y <-    pu[1]; x <-  0.5; hadj <- 0.5; vadj <- 1  } else
  if(titlepos==2) {y <- mean(pu); x <- -0.2; hadj <- 1  ; vadj <- 0.5} else
  if(titlepos==3) {y <-    pu[2]; x <-  0.5; hadj <- 0.5; vadj <- -0.2} else
  if(titlepos==4) {y <- mean(pu); x <-  1.2; hadj <- 0  ; vadj <- 0.5} else
  if(titlepos==5) {y <- mean(pu); x <-  0.5; hadj <- 0.5; vadj <- 0.5} else
  stop("Wrong titlepos. Must be integer between 1 and 5.")
  # actually write title:
  text(x=x, y=y, labels=title, adj=c(hadj, vadj), xpd=TRUE, ...)
    # kernel density:
  if(is.list(density) | is.null(density) | isTRUE(density) )
    {
    dp <- do.call(stats::density, args=owa(list(x=z), density))
    lines(y=dp$x, x=dp$y/max(dp$y))
    }
  } # end vertical -------------------------------------------------------------
  })
}
