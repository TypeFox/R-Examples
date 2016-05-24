## Copyright (C) 1995, 1996, 1997  Friedrich Leisch
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
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
## @deftypefn {Function File} {} fractdiff (@var{x}, @var{d})
## Compute the fractional differences @math{(1-L)^d x} where @math{L}
## denotes the lag-operator and @math{d} is greater than -1.
## @} # deftypefn

## Author: FL <Friedrich.Leisch@ci.tuwien.ac.at>
## Description: Compute fractional differences

fractdiff  <- function(x, d)  {

  N <- 100

  if (length(x) < 2)
    stop("fractdiff: x must be a vector")

  if (length(d) > 1)
    stop("fractdiff: d must be a scalar")

  if (d >= 1)
    for (k in 1:d)
      x <- x[-1] - x[-length(x)]

  if (d > -1) {

    ## Matlab rem returns negative output for negative input...
    sn <- d < 0
    d <- d %% 1
    if(sn) d <- d - 1 

    if (d != 0) {
      n <- 0:N
      w <- Re(gamma(-d+n) / gamma(-d) / gamma(n+1))
      retval <- fftfilt(w, x)
      retval <- retval[seq_along(x)]
    } else {
      retval <- x
    }

  } else {
    stop("fractdiff: d must be > -1")
  }
  retval
}
