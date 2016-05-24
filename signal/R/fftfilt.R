## Copyright (C) 1996, 1997 John W. Eaton
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
## @deftypefn {Function File} {} fftfilt (@var{b}, @var{x}, @var{n})
##
## With two arguments, @code{fftfilt} filters @var{x} with the FIR filter
## @var{b} using the FFT.
##
## Given the optional third argument, @var{n}, @code{fftfilt} uses the
## overlap-add method to filter @var{x} with @var{b} using an N-point FFT.
##
## If @var{x} is a matrix, filter each column of the matrix.
## @} # deftypefn

## Author: Kurt Hornik <Kurt.Hornik@ci.tuwien.ac.at>
## Created: 3 September 1994
## Adapted-By: jwe

   
FftFilter <- function(b, n) {
  res = list(b = b, n = n)
  class(res) = "FftFilter"
  res
}
 
filter.FftFilter <- function(filt, x, ...) 
  fftfilt(filt$b, x, filt$n)
  
fftfilt  <- function(b, x, n = NULL) {
 
  ## If n is not specified explicitly, we do not use the overlap-add
  ## method at all because loops are really slow.  Otherwise, we only
  ## ensure that the number of points in the FFT is the smallest power
  ## of two larger than N and length(b).  This could result in length
  ## one blocks, but if the user knows better ...
  N = n
  
  l_x = length(x)
  l_b = length(b)
 
  if (is.null(n)) {
    ## Use FFT with the smallest power of 2 which is >= length (x) +
    ## length (b) - 1 as number of points ...
    N = 2 ^ (ceiling(log(l_x + l_b - 1) / log(2)))
    B = fft(postpad(b, N))
    y = ifft(fft(postpad(x, N)) * B)
  } else {
    ## Use overlap-add method ...
    if (length(n) > 1)
      stop("fftfilt: n has to be a scalar")
    N = 2 ^ (ceiling(log(max(N, l_b)) / log(2)))
    L = N - l_b + 1
    B = fft(postpad(b, N))
#    B = B[,rep(1., l_x)]
    R = ceiling(l_x / L)
    y = numeric(l_x)
    for (r  in  1:R) {
      lo = (r - 1) * L + 1
      hi = min(r * L, l_x)
      tmp = numeric(0)
      tmp[1:(hi-lo+1)] = x[lo:hi]
      tmp = ifft(fft(postpad(tmp, N)) * B)
      hi  = min(lo+N-1, l_x)
      y[lo:hi] = y[lo:hi] + tmp[1:(hi-lo+1)]
    }
  }
  y = y[1:l_x]
 
  ## Final cleanups: if both x and b are real respectively integer, y
  ## should also be
 
  if (is.numeric(b) && is.numeric(x)) 
    y = Re(y)
  if (!any(as.logical(b - round(b)))) {
    idx = !any(as.logical(x - round(x)))
    y[idx] = round(y[idx])
  } 
  y
}
