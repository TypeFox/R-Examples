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

## usage: y = interp(x, q [, n [, Wc]])
##
## Upsample the signal x by a factor of q, using an order 2*q*n+1 FIR
## filter. Note that q must be an integer for this rate change method.
## n defaults to 4 and Wc defaults to 0.5.
##
## Example
##                                          # Generate a signal.
##    t=0:0.01:2; x=chirp(t,2,.5,10,'quadratic')+sin(2*pi*t*0.4)
##    y = interp(x(1:4:length(x)),4,4,1);   # interpolate a sub-sample
##    stem(t(1:121)*1000,x(1:121),"-g;Original;"); hold on
##    stem(t(1:121)*1000,y(1:121),"-r;Interpolated;")
##    stem(t(1:4:121)*1000,x(1:4:121),"-b;Subsampled;"); hold off
##
## See also: decimate, resample

interp <- function(x, q, n = 4, Wc = 0.5)  {

  if (q != round(q))
    stop("interp only works with integer q.")

  y = numeric(length(x)*q + q*n + 1)
  y[seq(1, length(x)*q, by = q)] = x
  b = fir1(2*q*n+1, Wc/q)
  y = q*fftfilt(b, y)
  y = y[-(1:(q*n+1))]  # adjust for zero filter delay
  y
}

#!demo
#! ## Generate a signal.
###t = seq(0, 2, by = 0.01)
###x = chirp(t,2,.5,10,'quadratic') + sin(2*pi*t*0.4)
###y = interp(x[seq(1, length(x), by = 4)],4,4,1)   # interpolate a sub-sample
###plot(t, x, type = "l")
###idx = seq(1,length(t),by = 4)
###lines(t, y[1:length(t)], col = "blue")
###points(t[idx], y[idx], col = "blue", pch = 19)
#! % graph shows interpolated signal following through the
#! % sample points of the original signal.
