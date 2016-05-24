#' Read data, so it can be used for high-risk zone methodology. 
#'
#' If xwin or ywin is NULL, the observation window will be a rectangular bounding box.
#' Vertices must be listed anticlockwise; no vertex should be repeated.
#' Only needed for data that is not already of class ppp.
#'
#' @param xppp  Vector of x coordinates of data points
#' @param yppp  Vector of y coordinates of data points
#' @param xwin  Vector of x coordinates of the vertices of a polygon circumscribing the observation window
#' @param ywin  Vector of y coordinates of the vertices of a polygon circumscribing the observation window
#' @param unitname Optional. Name of unit of length. Either a single character string, or a vector of two 
#' character strings giving the singular and plural forms, respectively.
#' @export  
#' @return An object of class "ppp" describing a point pattern in the two-dimensional plane.
#' @seealso \code{\link[spatstat]{ppp}}, \code{\link[spatstat]{bounding.box.xy}}, \code{\link[spatstat]{owin}}
#' @examples 
#' data(craterA)
#' windowA <- data.frame(x = craterA$window$bdry[[1]]$x, y = craterA$window$bdry[[1]]$y)
#' patternA <- data.frame(x = craterA$x, y = craterA$y)
#' str(patternA)
#' str(windowA)
#' crater <- read_pppdata(xppp = patternA$x, yppp = patternA$y, 
#'                        xwin = windowA$x, ywin = windowA$y)
#' crater

#- fr?her: readDataPPP
read_pppdata <- function(xppp, yppp, xwin=NULL, ywin=NULL, unitname=NULL) {
  
  #check if input arguments have correct values
  stopifnot(is.vector(xppp), is.vector(yppp))
  if ( !is.vector(xwin) & !is.null(xwin) ) stop("xwin must be a vector")
  if ( !is.vector(ywin) & !is.null(ywin) ) stop("ywin must be a vector")
  if ( is.null(xwin) | is.null(ywin) ) warning("since the coordinates for the window (xwin, ywin) are not given,
                                             it has to be calculated.")
  
  
  
  if (is.null(xwin) | is.null(ywin)){
    w <- bounding.box.xy(xppp, yppp)
    w$xrange <- w$xrange + c(-1,1)*0.1*(w$xrange[2] - w$xrange[1])
    w$yrange <- w$yrange + c(-1,1)*0.1*(w$yrange[2] - w$yrange[1])
    unitname <- spatstat::as.units(unitname)
    w$units <- unitname
  } else{
    win <- bounding.box.xy(xwin, ywin)
    w <- owin(win$xrange, win$yrange, poly=list(x=xwin, y=ywin), unitname=unitname)
  }
 
  
  result <- ppp(xppp, yppp, window=w)
  
  return(result)
}
