#' Points colored relative to third dimension
#' 
#' Draw colored points for 3D-data in a 2D-plane. Color is relative to third
#' dimension, by different classification methods. Can take 3 vectors or, as in
#' \code{\link{image}}, 2 vectors and a matrix for z.
#' 
#' @return Invisible list of values that can be passed to colPointsLegend or colPointsHist.
#' @note Rstudio scales graphics really badly, so don't expect the right legend width out of the box if you use Rstudio! 
#'      Exporting via \code{png("myplot.png", 600,400); colPoints(x,y,z); dev.off()} usually works much better
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2011-2014. I'd be interested in hearing what you used the function for.
#' @seealso \code{\link{classify}}, \code{\link{colPointsLegend}}, \code{\link{colPointsHist}}
#' @references \url{http://uxblog.idvsolutions.com/2011/10/telling-truth.html},
#'             \url{http://www.theusrus.de/blog/the-good-the-bad-22012/}
#' @keywords aplot hplot color
#' @export
#' @examples
#' 
#' i <- c( 22,  40,  48,  60,  80,  70,  70,  63,  55,  48,  45,  40,  30,  32)
#' j <- c(  5,  10,  15,  20,  12,  30,  45,  40,  30,  36,  56,  33,  45,  23)
#' k <- c(175, 168, 163, 132, 120, 117, 110, 130, 131, 160, 105, 174, 190, 183)
#' 
#' # basic usage:
#' colPoints(i,j,k, cex=1.5, pch="+", add=FALSE)
#' 
#' # with custom Range (only for method equalinterval):
#' colPoints(i,j,k, cex=1.5, pch="+", add=FALSE, Range=c(150, 180))
#' # can be used to allow comparison between several plots
#' # points outside the range are plotted with col2
#' 
#' # with custom colors:
#' mycols <- colorRampPalette(c("blue","yellow","red"))(50)
#' colPoints(i,j,k, cex=1.5, pch="+", add=FALSE, col=mycols)
#' 
#' # With legend title:
#' colPoints(i,j,k, cex=2, pch="+", add=FALSE,
#'          legargs=list(density=FALSE, title="Elevation [m above NN.]"))
#' ?colPointsLegend # to see which arguments can be set via legargs
#' 
#' # with lines (nint to change number of linear interpolation points):
#' colPoints(i,j,k, cex=1.5, add=FALSE, lines=TRUE, nint=10, lwd=2)
#' # With NAs separating lines:
#' tfile <- system.file("extdata/rivers.txt", package="berryFunctions")
#' rivers <- read.table(tfile, header=TRUE, dec=",")
#' colPoints(x,y,n, data=rivers, add=FALSE, lines=TRUE)
#' 
#' # different classification methods:
#' set.seed(007) ;  rx <- rnorm(30) ; ry <- rnorm(30) ; rz <- rnorm(30)*100
#' # sd: normal distribution
#' mycols <- colorRampPalette(c("blue","yellow", "red"))
#' colPoints(rx,ry,rz, add=FALSE, col=mycols(5), method="s",
#'           legargs=list(horiz=FALSE, x1=70, x2=95))
#' colPoints(rx,ry,rz, add=FALSE, col=mycols(6), method="s", sdlab=2,
#'           legargs=list(horiz=FALSE, labelpos=5, lines=FALSE, title=""))
#' # quantiles: each color is equally often used
#' colPoints(rx,ry,rz, add=FALSE, method="q",
#'           legargs=list(mar=c(0,5,3,0), bg="transparent") )
#' text(rx,ry,round(rz), col=8)
#' # logSpaced for rightly skewed data:
#' set.seed(41); rz2 <- rbeta(30, 1,7)*100
#' colPoints(rx,ry,rz2, add=FALSE, method="l", breaks=c(20,1.1708), col=mycols(20),
#'           legargs=list(mar=c(0,5,3,0), bg="transparent") )
#' colPoints(rx,ry,rz2, add=FALSE, method="q", breaks=0:20/20, col=mycols(20),
#'           legargs=list(mar=c(0,5,3,0), at=pretty2(rz2), labels=pretty2(rz2),
#'                        bg="transparent") )
#' 
#' # With histogram:
#' colPoints(i,j,k, add=FALSE, hist=TRUE)
#' colPoints(i,j,k, cex=3.5, lwd=3, pch=1, histargs=list(bg=5, breaks=5), add=FALSE)
#' colPoints(rx,ry,rz, cex=3.5, lwd=3, pch=1, add=FALSE, legend=FALSE,
#'    histargs=list(mar=c(0,0,0,0), x1=50,y1=99, x2=100,y2=80, yaxt="n"))
#' 
#' # use classify separately:
#' text(rx,ry,round(rz), col=mycols(100)[classify(rz)$index], cex=0.7)
#' 
#' # histogram in lower panel:
#' layout(matrix(1:2), heights=c(8,4) )
#' colPoints(i,j,k, add=FALSE, legargs=list(y2=80))
#' colPointsHist(z=k, x1=10,y1=80, x2=100,y2=10)
#' layout(1)
#' 
#' 
#' # Customizing the legend :
#' colPoints(i,j,k, legend=FALSE, add=FALSE)
#' colPointsLegend(x1=20,y1=50, x2=95,y2=40, z=k, labelpos=5, atminmax=TRUE, bg=7)
#' colPointsLegend(x1=50,y1=28, x2=90,y2=18, z=k, Range=c(80, 200), nbins=12, font=3)
#' colPointsLegend(x1=10,y1=15, x2=40,y2= 5, z=k, labelpos=5, lines=FALSE, title="")
#' 
#' colPointsLegend(z=k, horizontal=FALSE)
#' colPointsLegend(x1=1, y1=90, z=k, horizontal=FALSE, labelpos=4, cex=1.2)
#' colPointsLegend(x1=23,y1=95, z=k, horizontal=FALSE, labelpos=5, cex=0.8,
#'   dens=FALSE, title="", at=c(130,150,170), labels=c("y","rr","Be"), lines=FALSE)
#' # For method other than colPoints' default, it is easiest to include these
#' # options as a list in legargs, but you can also use the invisible output
#' # from colPoints for later calls to colPointsLegend
#' 
#' 
#' # colPoints with matrix:
#' colPoints(z=volcano, add=FALSE)
#' # image and contour by default transpose the matrix! This is really in the data
#' colPointsHist(z=volcano)
#' 
#' # highlight local character of points on a regular grid normally drawn with image:
#' # library(datasets), normally already loaded in newer R versions.
#' z <- t(volcano)  ;  x <- 1:ncol(z)  ;  y <- 1:nrow(z)
#' colPoints(x,y,z, add=FALSE)  # takes matrix for z
#' contour(x,y,t(z), add=TRUE)
#' 
#' # image only takes a regular matrix, but not scatterpoints...
#' image(x,y,t(z), col=rev(rainbow(100, start=0, end=.7)))
#' 
#' # add single newly measured points to image (fictional data):
#' mx <- c( 22,  40,  80,  45,  60,  63,  30,  70)
#' my <- c(  5,  33,  12,  56,  20,  40,  45,  45)
#' mz <- c(135, 155, 120, 105, 140, 130, 190, 110)
#' colPoints(mx,my,mz, cex=5, pch="*", Range=c(94, 195), col2=NA, legend=FALSE)
#' points(mx,my, cex=4)
#' text(mx,my,mz, adj=-0.5, font=2)
#' 
#' @param x,y Vectors with coordinates of the points to be drawn
#' @param z z values beloning to coordinates. Vector or matrix
#' @param data Optional: data.frame with the column names as given by x,y and z.
#' @param Range Ends of color bar for method=equalinterval. DEFAULT: range(z, finite=TRUE)
#' @param method Classification method (partial matching is performed), see \code{\link{classify}}. DEFAULT: "equalinterval")
#' @param breaks Specification for method, see \code{\link{classify}}. DEFAULT: different defaults for each method
#' @param sdlab Type of label and breakpoints if \code{method=standarddeviation}, see \code{\link{classify}}. DEFAULT: 1
#' @param col Vector of colors to be used. DEFAULT: \code{\link{rainbow}} from blue (lowest) to red (highest value in Range)
#' @param col2 Color for points where z is NA, lower or higher than Range. DEFAULT: c(NA, 1, 8)
#' @param legend Logical. Should a \code{\link{colPointsLegend}} be drawn? DEFAULT: TRUE
#' @param legargs List. Arguments passed to \code{\link{colPointsLegend}}. DEFAULT: NULL, with some defaults specified internally
#' @param hist Logical. Should a \code{\link{colPointsHist}} be drawn? DEFAULT: FALSE (TRUE if histargs are given)
#' @param histargs List. Arguments passed to \code{\link{colPointsHist}}. DEFAULT: NULL
#' @param add Logical. Should the points be added to current (existing!) plot? If FALSE, a new plot is started. DEFAULT: TRUE (It's called col\bold{Points}, after all)
#' @param lines Logical. Should lines be drawn underneath the points? (color of each \code{\link{segments}} is taken from starting point, last point is endpoint.) DEFAULT: FALSE
#' @param nint Numeric of length 1. Number of interpolation points between each coordinate if \code{lines=TRUE}. nint=1 means no interpolation. Values below 10 will smooth coordinates and miss the original points!. DEFAULT: 30
#' @param xlab x-axis label. DEFAULT: substitute as in plot
#' @param ylab y-axis label. DEFAULT: ditto
#' @param las Label Axis Style. Only used when add=FALSE. See \code{\link{par}}. DEFAULT: 1 (all labels horizontal)
#' @param pch Point CHaracter. See \code{\link{par}}. DEFAULT: 16
#' @param quiet Turn off warnings? DEFAULT: FALSE
#' @param \dots Further graphical arguments passed to plot, points and lines, eg cex, xlim (when add=F), mgp, main, sub, asp (when add=F), etc. Note: col does not work, as it is already another argument
#' 
colPoints <- function(
  x, y, # x,y: Vectors with coordinates of the points to be drawn
  z, # Vector or matrix with accompanying color defining height values
  data, # Optional: data.frame with the column names as given by x,y and z.
  Range=range(z, finite=TRUE), # Ends of color bar for method=equalinterval
  method="equalinterval", # type of binning or classification method (ways to get color class breakpoints)
  breaks, # specification for method
  sdlab=1, #
  col=seqPal(cl$nbins), # color palette. DEFAULT: 100 nuances from blue to red
  col2=c(NA, 1, 8), # color for z==NA and points not in the color range
  legend=TRUE, # Should a legend be drawn?
  legargs=NULL, # Arguments for colPointsLegend.
  hist=FALSE, # Should a legend be drawn?
  histargs=NULL, # Arguments for colPointsHist. FALSE to suppress drawing
  add=TRUE, # as in points. add to existing plot? add=F to draw new plot
  lines=FALSE, #  Logical. Should lines be drawn underneath the points?
  nint=30, # Numeric of length 1. Number of interpolation points between each coordinate if lines=TRUE.
  xlab=substitute(x), # axis labels
  ylab=substitute(y),
  las=1, # LabelAxisStyle: all labels horizontally (only relevant when add=FALSE)
  pch=16, # PointCHaracter, see ?par
  quiet=FALSE, # Turn off warnings?
  ...) # further arguments passed to plot, points and lines, eg cex, xlim (bei add=F), mgp, main, sub, asp (when add=F), etc. NOT col
{
xlab <- xlab ;  ylab <- ylab # defaults need to be set before x and y are evaluated
# error checking:
if(length(nint)>1) if(!quiet) warning("Only the first value of 'nint' is used.")
nint <- nint[1]
if(nint<1) stop("nint must be >= 1.")
col2 <- rep(col2, length.out=3) # in case only one, two or >3 values are given.
# Partial matching of method:
PossibleValues <- c("equalinterval", "quantile", "logspaced", "standarddeviation", "usergiven")
method <- PossibleValues[pmatch(tolower(method),  PossibleValues)]
if(is.na(method)) stop("method can only be equalinterval, quantile, logspaced, standarddeviation, or usergiven (but the name can be abbreviated).")
#
# vector vs matrix and dimension check: ----------------------------------------
# a) argument data is given
if(!missing(data)) # get x, y and z from data.frame
   {
   x <- data[ , deparse(substitute(x))]  
   y <- data[ , deparse(substitute(y))]  
   z <- data[ , deparse(substitute(z))] 
   } # now continue with case b
# error checking:
if(diff(range(z, finite=TRUE))==0) if(!quiet) warning("All z-values are equal.")
# b) Regular case: z ist a vector
if(is.vector(z))
   {
   if(!(length(x)==length(y) & length(x)==length(z)))
      stop("Vectors x,y,z are not all of the same length!")
   x <- x ;   y <- y ;   z <- z
   } else
# c) z is a matrix: class(z) = matrix, data.frame, array (2D) - as in image, persp
   {
   if(missing(x)) {x <- 1:ncol(z) ; xlab <- "x" }
   if(missing(y)) {y <- nrow(z):1 ; ylab <- "y" }
   if(!(length(x)==ncol(z) & length(y)==nrow(z)))
     stop("Dimension of z (ncol*nrow) is not length(x) * length(y)!")
   x <- rep(x, each=nrow(z));  y <- rep(y, ncol(z));  z <- as.vector(z)
   }
#
# CLASSIFICATION # -------------------------------------------------------------
if(method=="equalinterval") if(!missing(col)) breaks <- length(col)
#
cl <- classify(x=z, method=method, breaks=breaks, sdlab=sdlab, Range=Range, quiet=quiet)
# error check:
if(length(col) != cl$nbins) stop("Number of colors is not equal to number of classes.")
#
# ACTUAL PLOTTING --------------------------------------------------------------
if(!add) plot(x, y, col=NA, pch=pch, xlab=xlab, ylab=ylab, las=las, ...)
# Plot lines if wanted:
if(lines)
  {# linear interpolation between coordinates (smoother line colors):
  np <- length(x)*nint-nint+1 # replacing NA necessary if NAs are at start or end
  x2 <- approx(replace(x, is.na(x), median(x, na.rm=TRUE)), n=np)$y
  y2 <- approx(replace(y, is.na(y), median(y, na.rm=TRUE)), n=np)$y
  z2 <- approx(replace(z, is.na(z), median(z, na.rm=TRUE)), n=np)$y
  # classify interpolated values:
  cl2 <- classify(x=z2, method=method, breaks=breaks, sdlab=sdlab, Range=Range, quiet=quiet)
  # Where are NAs in the vectors?
  wNA <- is.na(x) | is.na(y) | is.na(z)
  # change single values (surrounded by NA) to NA:
  for(i in 2:(length(wNA)-1))
      if(isTRUE(wNA[i-1]) & isTRUE(wNA[i+1]))  wNA[i] <- TRUE
  # change interpolated values to NA where corresponding values are NA:
  for(i in which(wNA))
      cl2$index[pmax((i-2)*nint+1, 1) : pmin(i*nint, np)] <- NA
  # Actually draw segments:
  segments(x0=x2[-length(x2)],  y0=y2[-length(y2)],  x1=x2[-1],  y1=y2[-1],
           col=col[cl2$index], ...)
  }
points(x[is.na(z)], y[is.na(z)], col=col2[1], pch=pch, ...)
points(x, y, col=c(col, col2[2:3])[cl$index], pch=pch, ...)
#
# add legend:
if(legend)
  {
  legdefs <- list(z=z, at=cl$at, labels=cl$labels, bb=cl$bb, nbins=cl$nbins,
      colors=col, plottriangle=any(na.omit(cl$index>cl$nbins)), tricol=col2[2:3])
  do.call(colPointsLegend, args=owa(legdefs, legargs))
  }
#
# add histogramm:
if(hist | !missing(histargs))
  {
  histdefs <- list(z=z, at=cl$at, labels=cl$labels, bb=cl$bb, nbins=cl$nbins, colors=col)
  do.call(colPointsHist, args=owa(histdefs, histargs))
  }
return(invisible(cl))
} # Function end
