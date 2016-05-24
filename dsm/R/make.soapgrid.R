#' Create a knot grid for the internal part of a soap film smoother.
#'
#' This routine simply creates a grid of knots (in the correct format) to be
#' used as in the "internal" part of the soap film smoother
#'
#' @param bnd list with elements \code{x} and \code{y} which give the locations
#'        of the boundary vertices. The first and last elements should be the
#'        same.
#' @param n.grid either one number giving the number of points along the 
#'        \code{x} and \code{y} axes that should be used to create the grid, or
#'        a vector giving the number in the \code{x} direction, then \code{y}
#'        direction.
#'
#' @return a list with elements \code{x} and \code{y}, containing the knot 
#'         locations.
#'
#' @author David L Miller
#' @export
#'
make.soapgrid<-function(bnd,n.grid){
  # set the grid size, if the input is a 2-vec then it is m and n
  if(length(n.grid)==2){
     m<-n.grid[1]
     n<-n.grid[2]
  }else{
     m<-n<-n.grid
  }

  # min and max values of the boundary (but not on the boundary)
  xmin<-min(bnd$x,na.rm=TRUE)
  ymin<-min(bnd$y,na.rm=TRUE)
  xmax<-max(bnd$x,na.rm=TRUE)
  ymax<-max(bnd$y,na.rm=TRUE)

  # create the grid
  ng<-expand.grid(x=seq(xmin,xmax,length=m),
                  y=seq(ymin,ymax,length=n))

  # which bits of the grid were inSide?
  x <- ng$x; y <- ng$y
  onoff<-inSide(as.list(bnd),x,y)

  # remove the outside ones
  ng<-ng[onoff,]

  return(ng)
}
