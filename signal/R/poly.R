## Copyright (C) 1996, 1997 John W. Eaton
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it andsolve(,)or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, write to the Free
## Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.

## -*- texinfo -*-
## @deftypefn {Function File} {} poly (@var{a})
## If @var{a} is a square @math{N}-by-@math{N} matrix, @code{poly (@var{a})}
## is the row vector of the coefficients of @code{det (z * eye (N) - a)},
## the characteristic polynomial of @var{a}.  If @var{x} is a vector,
## @code{poly (@var{x})} is a vector of coefficients of the polynomial
## whose roots are the elements of @var{x}.
## @} # deftypefn

## Author: KH <Kurt.Hornik@neuro.tuwien.ac.at>
## Created: 24 December 1993
## Adapted-By: jwe

poly <- function(x) {
  n = NROW(x)
  m = NCOL(x)
  if (is.null(x) || length(x) == 0)
    return(1)
  if (m == 1) {
    v = x
  } else if (m == n) {
    v = eigen(x)$values
  } else {
    stop("poly(x), where x is a vector or a square matrix")
  }
  
  y = numeric(n+1)
  y[1] = 1
  for (j in 1:n) 
    y[2:(j+1)] = y[2:(j+1)] - v[j] * y[1:j]

  if (all(Im(x) == 0))
    y = Re(y)

  y
}
