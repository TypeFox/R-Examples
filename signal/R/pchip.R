## Copyright (C) 2001,2002  Kai Habel
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  
## -*- texinfo -*-
## @deftypefn {Function File} {@var{pp} = } pchip (@var{x}, @var{y})
## @deftypefnx {Function File} {@var{yi} = } pchip (@var{x}, @var{y}, @var{xi})
## piecewise cubic hermite interpolating polynom.
## @var{x} must be a strictly monotonic vector (either increasing or decreasing).
## @var{y} is a vector of same length as @var{x},
## or a matrix where the number of columns must match the length
## of @var{x}. In this case the interpolating polynoms are calculated
## for each column. 
## In contrast to spline, pchip preserves the monotonicity of (@var{x},@var{y}).
## 
## @seealso{ppval, spline, csape}
## @} # deftypefn
  
## Author:  Kai Habel <kai.habel@gmx.de>
## Date: 9. mar 2001
##
## S_k = a_k + b_k*x + c_k*x^2 + d_k*x^3; (spline polynom)
##
## 4 conditions:
## S_k(x_k) = y_k
## S_k(x_k+1) = y_k+1
## S_k'(x_k) = y_k'
## S_k'(x_k+1) = y_k+1'
  
pchip  <- function(x, y, xi = NULL)  {
  
  pchip_deriv <- function(x, y) {
    if (NROW(x) != NROW(y))
      stop("number of rows for x and y must match")
    # Octave does this for each column of y; we just assume a vector.
    z <- .Fortran("dpchim",
                  as.integer(NROW(x)),
                  as.double(x),
                  as.double(y),
                  output = as.double(y),
                  as.integer(1),
                  ierr = as.integer(1),
                  PACKAGE = "signal")
    if (z$ierr >= 0)
      return(z$output)
    else
      stop("error")
  }
  
  n = length(x)
  
  ry = NROW(y)
  cy = NCOL(y)
  
  h = diff(x)
  if (all(h < 0)) {
    x = rev(x)
    h = diff(x)
    y = as.matrix(y[, cy:1])
  } else if (any(h<=0)) {
    stop('x must be strictly monotonic')
  }
  
  if (ry != n)
    stop("size of x and y must match")
  
  if (cy > 1)
    h = diff(x) %x% rep.int(1, cy) # kron multiplication
  
  ## not used: ## dy = diff(y) / h
  
  a = y
  b = pchip_deriv(x,y)
  c = - (b[2:n] + 2 * b[1:(n-1)]) / h + 3 * diff(a) / h ^ 2
  d = (b[1:(n-1)] + b[2:n]) / h^2 - 2 * diff(a) / h^3
  
  d = d[1:(n-1)]; c = c[1:(n-1)]
  b = b[1:(n-1)]; a = a[1:(n-1)]
  coeffs = cbind(d, c, b, a)
  pp = mkpp(x, coeffs)
  
  if (is.null(xi))
    ret = pp
  else
    ret = ppval(pp,xi)
  ret
} 
  
#xf = seq(0,10, length=500); yf = sin(2*pi*xf/5)
#xp = 0:10; yp = sin(2*pi*xp/5)
#xp = c(0:1,3:10); yp = sin(2*pi*xp/5)
#lin  = pchip(xp, yp, xf)
#lin  = interp1(xp, yp, xf, 'pchip')
#plot(xp, yp)
#lines(xf, lin, col = "red")
  
  
## Copyright (C) 2000 Paul Kienzle
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  
## usage: pp = mkpp(x, P)
## Construct a piece-wise polynomial structure from sample points x and
## coefficients P.  The ith row of P, P(i,:), contains the coefficients 
## for the polynomial over the ith interval, ordered from highest to 
## lowest. There must be one row for each interval in x, so 
## rows(P) == length(x)-1.  
##
## You can concatenate multiple polynomials of the same order over the 
## same set of intervals using P = [ P1 ; P2 ; ... ; Pd ].  In this case,
## rows(P) == d*(length(x)-1).
##
## mkpp(x, P, d) is provided for compatibility, but if d is not specified
## it will be computed as round(rows(P)/(length(x)-1)) instead of defaulting
## to 1.
  
mkpp <- function(x, P, d = round(NROW(P)/pp$n)) {
  pp = list()
  pp$x = x
  pp$P = P
  pp$n = length(x) - 1
  pp$k = NCOL(P)
  pp$d = d
  if (pp$n*d != NROW(P))
    stop("mkpp: num intervals in x doesn't match num polynomials in P")
  class(pp) = "pp"
  pp
} 
  
##%!demo # linear interpolation
#x = seq(0,pi,length=5) 
#t= cbind(sin(x), cos(x))
#m = diff(t) / (x[2]-x[1])
#b = t[1:4,]
#pp = mkpp(x, cbind(as.vector(m),as.vector(b)))
#ppval(pp,x)
##%! xi=linspace(0,pi,50)
##%! plot(x,t)
##   lines(xi,ppval(pp,x))

## Copyright (C) 2000 Paul Kienzle
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## usage: yi = ppval(pp, xi)
## Evaluate piece-wise polynomial pp and points xi.  If pp.d > 1,
## the returned yi will be a matrix of size rows(xi) by pp.d, or
## its transpose if xi is a row vector.

ppval = function(pp, xi) {
  if (is.null(xi)) {
    yi = NULL
  } else {
    lookup = function(x, xi) approx(x, 1:length(x), xi, method="constant", yleft=0, yright=length(x))$y
    idx = lookup(pp$x[2:pp$n], xi) + 1
    dx = as.matrix(xi - pp$x[idx])
    dx = dx[,rep(1,pp$d)]
    c = matrix(pp$P, pp$n, pp$d)
    yi = c[idx, ]
    if (pp$k > 1)
      for (i  in 2:pp$k) {
        c = matrix(pp$P[, i], pp$n, pp$d)
        yi = yi * dx + c[idx, ]
      }
  }
  yi
}
