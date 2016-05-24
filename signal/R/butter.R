## Copyright (C) 1999 Paul Kienzle
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

## Generate a butterworth filter.
## Default is a discrete space (Z) filter.
## 
## [b,a] = butter(n, Wc)
##    low pass filter with cutoff pi*Wc radians
##
## [b,a] = butter(n, Wc, 'high')
##    high pass filter with cutoff pi*Wc radians
##
## [b,a] = butter(n, [Wl, Wh])
##    band pass filter with edges pi*Wl and pi*Wh radians
##
## [b,a] = butter(n, [Wl, Wh], 'stop')
##    band reject filter with edges pi*Wl and pi*Wh radians
##
## [z,p,g] = butter(...)
##    return filter as zero-pole-gain rather than coefficients of the
##    numerator and denominator polynomials.
## 
## [...] = butter(...,'s')
##     return a Laplace space filter, W can be larger than 1.
## 
## [a,b,c,d] = butter(...)
##  return  state-space matrices 
##
## References: 
##
## Proakis & Manolakis (1992). Digital Signal Processing. New York:
## Macmillan Publishing Company.

## Author: Paul Kienzle <pkienzle@user.sf.net>
## Modified by: Doug Stewart <dastew@sympatico.ca> Feb, 2003

butter <- function(n, ...) UseMethod("butter")

butter.FilterOfOrder <- function(n, ...)
  butter(n$n, n$W, n$type, ...)

butter.default <- function(n, W, type = c("low", "high", "stop", "pass"), plane = c("z", "s"), ...) {

  type <- match.arg(type)
  plane <- match.arg(plane)

  ## interpret the input parameters
  if (!(length(n)==1 && n == round(n) && n > 0))
    stop("butter: filter order n must be a positive integer")

  stop <- type == "stop" || type == "high"
  digital <- plane == "z"

  if (length(W) != 1 && length(W) != 2)
    stop("butter: frequency must be given as w0 or c(w0, w1)")

  if (digital && !all(W >= 0 & W <= 1))
    stop("butter: critical frequencies must be in (0 1)")
  else if (!digital && !all(W >= 0))
    stop("butter: critical frequencies must be in (0 inf)")

  ## Prewarp to the band edges to s plane
  if (digital) {
    T <- 2       # sampling frequency of 2 Hz
    W <- 2 / T*tan(pi * W / T)
  }

  ## Generate splane poles for the prototype butterworth filter
  ## source: Kuc
  C <- 1 # default cutoff frequency
  pole <- C*exp(1i*pi*(2*1:n + n - 1) / (2*n))
  if (n %% 2 == 1)
    pole[(n+1) / 2] <- -1  # pure real value at exp(i*pi)
  zero <- numeric(0)
  gain <- C^n

  ZPG <- Zpg(zero = zero, pole = pole, gain = gain)

  ## s-plane frequency transform
  ZPG <- sftrans(ZPG, W = W, stop = stop)

  ## Use bilinear transform to convert poles to the z plane
  if (digital)
     ZPG <- bilinear(ZPG, T = T)

  as.Arma(ZPG)
} 
