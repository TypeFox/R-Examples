## Copyright (C) 2015  Szymon Sacher, Andrew Clausen The University of Edinburgh
##
## Excerpts adapted from F77 code Copyright (C) Paolo Costantini
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

## This program is based on results from paper below and on supplementary
## F77 code obtained from the authors. Questions should be sent to maintainer,
## Szymon Sacher at s1340144@sms.ed.ac.uk. Relevant pieces of original paper
## are cited in [].
##
## 1/ Costantini, P; Fontanella, F; 'Shape-Preserving Bivariate Interoplation',
## SIAM J Numer. Anal. Vol.27, No. 2, pp. 488-506, 1990

sps_fun <- function(x, y, z=NULL, der=0, der.x=der, der.y=der, grid = FALSE, ...)

{
  dim <- ifelse(is.null(z), 1, 2)

  spline <- sps_prep(x, y, z, ...)

  if (dim == 1)
    function(x) sps_eval(spline, x, der.x)
  else
    function(x,y) sps_eval(spline, x, der.x, y, der.y, grid)
}

sps_interpolate <- function(x, y, z=NULL, xt, yt=NULL, grid=TRUE, ...)

{
  stopifnot(is.null(z) == is.null(yt))

  fun <- sps_fun(x, y, z, grid=grid, ...)

  dim <- ifelse(is.null(z), 1, 2)

  if (dim == 1)
    fun(xt)
  else
    fun(xt, yt)
}

sps_prep <- function(
  x, y, z=NULL,
  fx = NA, fy = NA, fxy = NA,
  shape = c("monotonicity", "curvature"),
  shape.x = shape, shape.y = shape,
  max.deg = 50, smoothness = 1,
  tol = 0.0001)
{


  if (round(smoothness, 0) != smoothness || smoothness < 1)
    stop('"smoothness" must be positive integer')

  if (!all(c(shape, shape.x, shape.y) %in% c('monotonicity', 'curvature')))
    stop('shape, shape.x, and shape.y must only contain "monotonicity" and/or "curvature"')



  if(is.null(z)){

    stopifnot(length(x) == length(y))

    stopifnot(is.vector(x) && is.vector(y))
    spline <- sps_prep_u(
      x, y, fx, shape.x, max.deg, smoothness, tol)
  }

  else{
    stopifnot(length(x) == length(z) || all(dim(z) == c(length(x),length(y))))

    if(length(x) == length(z)) {
      gridded <- as.grid(x,y,z)
      x <- gridded$x
      y <- gridded$y
      z <- gridded$z
    }

    spline <- sps_prep_bi(
      x,y,z,fx,fy,fxy,shape.x,shape.y,max.deg,smoothness,tol)
  }
  spline
}

# Attempts to transform triples (x,y,z) to grid form so that x and y are vectors of
# coordinates and z is matrix of values no grid spanned by x and y.

as.grid <- function(x, y, z)
{
  x_grid <- sort(unique(x))
  y_grid <- sort(unique(y))

  N <- length(x_grid)
  M <- length(y_grid)

  gridded <- matrix(ncol=M, nrow=N)

  for(k in 1:length(x)){
    i <- which(x_grid == x[k])
    j <- which(y_grid == y[k])

    if(!is.na(gridded[i,j])) stop('Unable to transform to grid. Repeated data.')
    gridded[i,j] <- z[k]
  }

  if(any(is.na(gridded))) stop('Unable to transform to grid. Missing data.')

  list(x = x_grid, y=y_grid, z=gridded)
}

sps_eval<- function(spline, x, der.x=NULL, y=NULL, der.y=NULL, grid=FALSE)

{
  if (spline$dim == 1)
    sps_eval_u(x,spline,der.x)
  else if (spline$dim == 2)
    sps_eval_bi(x,y,spline,der.x,der.y, grid)
  else
    stop('Invalid spline$dim')

}




