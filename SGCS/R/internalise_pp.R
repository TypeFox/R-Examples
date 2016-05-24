#' Convert input to SGCS internal point pattern data
#' 
#' @param x some form of point pattern data.
#' 
#' @details
#' Understands 'ppp', 'pp3', 'matrix' and list(coordinate_matrix, bbox).
#' 
#' Data.frame marks are not well understood, so convert the mark data to
#' suit the analysis before use. 
#'  
#' Mainly for internal use.
#' 
#' @export
internalise_pp <- function(x) {
  dim <- 2
  ### 2D ppp:
  if(class(x)=="ppp"){
    coord <- as.matrix(coords(x))
    bbox <- with(x$window, cbind(xrange, yrange))
    pp <- list(coord=coord, n=nrow(coord), dim=ncol(coord), bbox=bbox, owin=x$window)
    if(!is.null(x$marks)) {
      if(x$markformat != "vector") stop("data.frame marks not yet supported.")
      if(is.factor(x$marks)) pp$type <- as.integer(x$marks)
      else pp$mass <- as.numeric(x$marks)
    }
    # add window area
    pp$area <- area(x$window) 
  }
  #### spatstat 3d pattern
  else if(is.pp3(x)){ 
    pp <- internalise_pp(list(coord=as.matrix(coords(x)), 
                                      bbox = with(x$domain, cbind(xrange, yrange, zrange)))
                              )
                         
  }
  #### if list of coordinates and bounding box -format
  else if(is.list(x)){
    if(is.null(x$coord) & !is.null(x$x)) if(ncol(x$x)>1) x$coord <- x$x
    if(is.null(x$coord)) stop("Can not interpret input pattern.")
    bbox <- x$bbox
    coord <- x$coord
    n <- nrow(coord)
    dim <- ncol(coord)
    pp <- list(coord = coord, bbox=bbox, n=n, dim=dim)
    # check for marks
    if(length(x$mark)==n){
      if(is.factor(x$mark)) pp$type <- as.integer(x$mark)
      else pp$mass <- as.numeric(x$marks)
    }
    if(length(x$mass)==n) pp$mass <- as.numeric(x$mass)
    if(length(x$type)==n) pp$type <- as.integer(x$type)
    pp$dim <- dim
    pp$n   <- nrow(pp$coord)
    pp$area <- prod( apply(pp$bbox, 2, diff) )#
  }
  else if(is.matrix(x) | is.data.frame(x)){ # matrix of coordinates given
    coord <- as.matrix(x)
    bbox <- bounding_box_xy(x)
    pp <- list(coord = coord, bbox=bbox, n=nrow(coord), dim=ncol(coord))
    pp$dim <- dim
    pp$n   <- nrow(pp$coord)
    pp$area <- prod( apply(pp$bbox, 2, diff) )
  }
  else stop("Can not interpret input pattern.")
  
  if(! pp$dim %in% c(2,3)) stop("Only 2- and 3- column matrices of coordinates supported.")
  ## some failsafes
  pp$bbox <- pp$bbox + 0.0
  pp$area <- pp$area + 0.0
  pp$mass <- pp$mass + 0.0
  ## done
  pp
}
