## Copyright (C) 2006 EPRI Solutions, Inc.
## by Tom Short, tshort@eprisolutions.com
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

remez <- function(n, f, a, w = rep(1.0, length(f) / 2),
                  ftype = c('bandpass', 'differentiator', 'hilbert'),
                  density = 16) {

  ftype = as.integer(factor(match.arg(ftype), c('bandpass', 'differentiator', 'hilbert')))

  if (n < 4)
    stop("Number of taps (n) must be an integer greater than 3")
  if (length(f) %% 2 == 1) 
    stop("'f' must be even")
  if (any(diff(f) < 0))
    stop("'f' must be nondecreasing")
  if (any(f < 0) || any(f > 1))
    stop("'f' must be in the range [0,1]")
  if (length(a) != length(f)) 
    stop("length(a) must equal length(f)")
  if (2*length(w) != length(f)) 
    stop("length(w) must be half of length(f)")
  if (density < 16)
    stop("density is too low, must be greater than or equal to 16")

  z <- .C("remez",
          coeff = as.double(rep(0, n + 1)),
          as.integer(n+1),
          as.integer(length(f) / 2),
          as.double(f / 2),
          as.double(a),
          as.double(w),
          as.integer(ftype),
          as.integer(density),
          PACKAGE = "signal")

  Ma(z$coeff)
}

#remez(15,c(0,0.3,0.4,1),c(1,1,0,0))
