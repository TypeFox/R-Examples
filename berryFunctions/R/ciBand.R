#' polygon confidence bands
#' 
#' \code{\link{polygon}} for confidence interval bands, can handle NA's well
#' 
#' @return None, currently. Used for drawing.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, July 2015
#' @seealso \code{\link{quantileBands}}, \code{\link{polygon}}, \code{\link{approx2}}
#' @keywords aplot hplot
#' @export
#' @examples
#' 
#' y1 <- c(1,3,4,2,1,4,6,8,7)
#' y2 <- c(5,6,5,6,9,8,8,9,10)
#' y3 <- c(4,4,5,4,4,6,7,8,9)
#' ciBand(yl=y1, yu=y2, ym=y3)
#' 
#' y1[6:7] <- NA
#' ciBand(yl=y1, yu=y2, ym=y3) # interpolation marked with stars if nastars=TRUE
#' ciBand(yl=y1, yu=y2, ym=y3, na="remove")
#' lines(y1, col=3, type="o")
#' lines(y2, col=3, type="o")
#' 
#' y2[1] <- NA
#' ciBand(yl=y1, yu=y2, ym=y3) # next observation carried backwards (NAs at begin)
#' # LOCF (last observation carried forwards if NAs at end)
#' # See ?approx2 for median/mean imputation in these cases
#' ciBand(yl=y1, yu=y2, ym=y3, na="remove")
#' y2[9] <- NA
#' ciBand(yl=y1, yu=y2, ym=y3)
#' ciBand(yl=y1, yu=y2, ym=y3, na="remove")  # NAs at both ends
#' y2[1] <- 5
#' ciBand(yl=y1, yu=y2, ym=y3)
#' ciBand(yl=y1, yu=y2, ym=y3, na="remove")  # NA only at end
#' 
#' # Actual usefull stuff:  sample size dependency of max and mean
#' ssdep_max <- function(n) quantile(  replicate(n=200, expr=max(rnorm(n)) )  )
#' ssdep_mean<- function(n) quantile(  replicate(n=200,expr=mean(rnorm(n)) )  )
#' x <- 1:100
#' res_max <- sapply(x, ssdep_max)
#' res_mean <- sapply(x, ssdep_mean)
#' ciBand(yl=res_max[2,], yu=res_max[4,], ym=res_max[3,], x=x, ylim=c(-0.5, 3))
#' ciBand(res_mean[2,], res_mean[4,], res_mean[3,], x=x, add=TRUE, colm="purple")
#' 
#' @param yu y values of upper confidence region boundary
#' @param yl y values of lower confidence region boundary
#' @param ym y values of median/mean line. Only added if this argument is given. DEFAULT: NULL
#' @param x x values (one ascending vector). DEFAULT: 1:length(yu)
#' @param na Method used at NA points. One of "interpolate" or "remove". DEFAULT: "interpolate"
#' @param nastars If na="interpolate", should stars be drawn at places that used to be NA? DEFAULT: TRUE
#' @param singlepoints If na="remove", add points for places surrounded by NAs? 
#'       can be a boolean (T/F) vector of length three for upper, lower, median. 
#'       Code to identify isolated points is taken from wq::plotTs. DEFAULT: TRUE
#' @param args List of arguments passed to \code{\link{points}} for the previous two arguments. DEFAULT: NULL
#' @param add Add to existing plot? If FALSE, plot is called before adding confidence interval. DEFAULT: FALSE
#' @param colm Color for median/mean line. DEFAULT: "green3"
#' @param colb Color of the confidence region band. DEFAULT: addAlpha(colm)
#' @param border \code{\link{polygon}} border. DEFAULT: NA
#' @param las LabelAxisStyle (axis labels turned upright, see \code{\link{par}}). DEFAULT: 1
#' @param ylim limits of plot. DEFAULT: range(yu,yl, finite=TRUE)
#' @param \dots Further arguments passed to \code{\link{plot}} - or maybe better polygon??
#' 
ciBand <- function(
yu,                  
yl,                  
ym=NULL,            
x=1:length(yu),      
na="interpolate",    
nastars=TRUE,        
singlepoints=TRUE,   
args=NULL,           
add=FALSE,          
colm="green3",      
colb=addAlpha(colm), 
border=NA,           
las=1,              
ylim=range(yu,yl, finite=TRUE), 
...   
)
{
# input checking:
# checks to control if inputs can be converted to vectors?
# ...
# recycle singlepoints.
singlepoints <- rep(singlepoints, length.out=3)
# Vector length checking:
nyl <- length(yl)
if( length(yu) != nyl) stop("Vectors yu and yl are not of the same length. (",
                         length(yu), " and ", nyl, " ).")
 if(!is.null(ym)) if( length(ym) != nyl ) stop("Vectors ym and yu/yl are not of the same length. (",
                         length(ym), " and ", nyl, " ).")
if( length( x) != nyl ) stop("Vectors x and yu/yl are not of the same length. (",
                         length( x), " and ", nyl, " ).")
if(!add) plot(x, yu, type="n", las=1, ylim=ylim, ...)
if(na=="interpolate")
  {
  # Interpolations:
  yui <- approx2(yu)
  yli <- approx2(yl)
  x <- approx2(x)
  if(!is.null(ym)) ymi <- approx2(ym)
  # Actual polygon:
  polygon(x=c(x, rev(x)), y=c(yui, rev(yli)), col=colb, border=border)
  # NA stars:
  if(nastars) do.call(points, args=owa(list(x=x[is.na(yu)], y=yui[is.na(yu)], pch=8), args, "x", "y"))
  if(nastars) do.call(points, args=owa(list(x=x[is.na(yl)], y=yli[is.na(yl)], pch=8), args, "x", "y"))
  # Draw median/mean line:
  if(!is.null(ym)) lines(x, ymi, col=colm)
  }
else if(na=="remove")
  {
  nna <- is.na(yu) | is.na(yl) | is.na(x)
  # break vectors up into sections of consecutive non-NA values:
  r <- rle(nna)
  streaks <- cumsum(r$lengths)
  nonna <- which(!r$values)
  sectionstarts <- c(1,streaks+1)[nonna]
  sectionends   <- streaks[nonna]
  #sectionstarts <- c(1, which(diff(nna)==-1)+1)
  #sectionends   <- c(which(diff(nna)==1), nyl)
  if(length(sectionstarts)!=length(sectionends)) stop("non-NA Sections not identified correctly.")
  # Actual polygons:
  for(i in 1:length(sectionstarts))
    {
    use <- sectionstarts[i]:sectionends[i]
    yur <- yu[use]
    ylr <- yl[use]
    xr  <-  x[use]
    polygon(x=c(xr, rev(xr)), y=c(yur, rev(ylr)), col=colb, border=border)
    }
  # Single points (isolated points, surrounded by NA)
  if(any(singlepoints))
    {
    iso <- function(x) # Code taken from wq::plotTs
      {
      x.forward <- c(NA, x[1:(length(x)-1)])
      x.back <- c(x[2:length(x)], NA)
      iso.pts <- is.na(x.forward) & is.na(x.back) & !is.na(x)
      #iso <- ifelse(iso.pts, x, NA)
      iso.pts
      }
    if(singlepoints[1]) do.call(points, args=owa(list(x=x[iso(yu)], y=yu[iso(yu)], 
                                             pch=20, col=colb), args, "x", "y"))
    if(singlepoints[2]) do.call(points, args=owa(list(x=x[iso(yl)], y=yl[iso(yl)], 
                                             pch=20, col=colb), args, "x", "y"))
    }
  # Draw median/mean line:
  if(!is.null(ym)) 
    {
    lines(x, ym, col=colm)
    if(singlepoints[3]) do.call(points, args=owa(list(x=x[iso(yl)], y=yl[iso(yl)], 
                                             pch=20, col=colm), args, "x", "y"))
    }
  } # End of na="remove"
else stop("na method ", na, " is not available (yet).")
}

