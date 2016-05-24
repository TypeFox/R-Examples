## Copyright (C) 1999, 2006, 2007, 2009 David M. Doolin
## bugs and limitations:  
##        Probably ought to be an optional check to make sure that
##        traversing the vertices doesn't make any sides cross 
##        (Is simple closed curve the technical definition of this?). 

##' Determines area of a polygon by triangle method.  The variables
##' \code{x} and \code{y} define the vertex pairs, and must therefore
##' have the same shape.  They can be either vectors or arrays.  If
##' they are arrays then the columns of \code{x} and \code{y} are
##' treated separately and an area returned for each.
##'
##' If the optional \code{dim} argument is given, then \code{polyarea}
##' works along this dimension of the arrays \code{x} and \code{y}.
##'
##' @title Determines area of a polygon by triangle method. 
##' @param x X coordinates of verticies.
##' @param y Y coordinates of verticies.
##' @param d Dimension of array to work along.
##' @return Area(s) of polygon(s).
##' @author David Sterratt based on the octave sources by David M. Doolin
##' @export
##' @import magic
##' @examples
##' x <- c(1, 1, 3, 3, 1)
##' y <- c(1, 3, 3, 1, 1)
##' polyarea(x, y)
##' polyarea(cbind(x, x), cbind(y, y)) ##  c(4, 4)
##' polyarea(cbind(x, x), cbind(y,y), 1) ##  c(4, 4)
##' polyarea(rbind(x, x), rbind(y,y), 2) ##  c(4, 4)
polyarea <- function(x, y, d=1) {
  if (is.vector(x) & is.vector(y)) {
    if (length(x) == length(y)) {
      a <- abs(sum(x*(shift(y, -1) - shift(y, 1))))/2
    } else {
      stop("x and y must have the same shape")
    }
  } else {
    if (is.array(x) & is.array(y)) {
      if (length(dim(x)) != 2) {
        stop("Arrays must have two dimensions.")
      }
      if (all(dim(x) == dim(y))) {
        v <- c(0, 0)
        v[d] <- 1
        a <- abs(apply(x*(ashift(y, -v) - ashift(y, v)), 3 - d, sum))/2
      } else {
        stop("x and y must have the same shape")
      }
    } else {
      stop("x and y must be of same type")
    }
  }
  return(a)
}

## %!shared x, y
## %! x = [1;1;3;3;1];
## %! y = [1;3;3;1;1];
## %!assert (polyarea(x,y), 4, eps)
## %!assert (polyarea([x,x],[y,y]), [4,4], eps)
## %!assert (polyarea([x,x],[y,y],1), [4,4], eps)
## %!assert (polyarea([x,x]',[y,y]',2), [4;4], eps)
