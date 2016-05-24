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

## Generate an Chebyshev type I filter with Rp dB of pass band ripple.
## 
## [b, a] = cheby1(n, Rp, Wc)
##    low pass filter with cutoff pi*Wc radians
##
## [b, a] = cheby1(n, Rp, Wc, 'high')
##    high pass filter with cutoff pi*Wc radians
##
## [b, a] = cheby1(n, Rp, [Wl, Wh])
##    band pass filter with edges pi*Wl and pi*Wh radians
##
## [b, a] = cheby1(n, Rp, [Wl, Wh], 'stop')
##    band reject filter with edges pi*Wl and pi*Wh radians
##
## [z, p, g] = cheby1(...)
##    return filter as zero-pole-gain rather than coefficients of the
##    numerator and denominator polynomials.
##
## [...] = cheby1(...,'s')
##     return a Laplace space filter, W can be larger than 1.
## 
## [a,b,c,d] = cheby1(...)
##  return  state-space matrices 
## 
## References: 
##
## Parks & Burrus (1987). Digital Filter Design. New York:
## John Wiley & Sons, Inc.

## Author: Paul Kienzle <pkienzle@user.sf.net>
## Modified: Doug Stewart Feb. 2003

cheby1 <- function(n, ...) UseMethod("cheby1")

cheby1.FilterOfOrder <- function(n, Rp = n$Rp, W = n$Wc, type = n$type, ...)
  cheby1(n$n, Rp, W, type, ...)

cheby1.default <- function(n, Rp, W, type = c("low", "high", "stop", "pass"), plane = c("z", "s"), ...) {

  type <- match.arg(type)
  plane <- match.arg(plane)

  ## interpret the input parameters
  if (!(length(n)==1 && n == round(n) && n > 0))
    stop("cheby1: filter order n must be a positive integer")

  stop <- type == "stop" || type == "high"
  digital <- plane == "z"

  if (length(W) != 1 && length(W) != 2)
    stop("cheby1: frequency must be given as w0 or c(w0, w1)")

  if (digital && !all(W >= 0 & W <= 1))
    stop("cheby1: critical frequencies must be in (0 1)")
  else if (!digital && !all(W >= 0))
    stop("cheby1: critical frequencies must be in (0 inf)")

  if (Rp < 0)
    stop("cheby1: passband ripple must be positive decibels")

  ## Prewarp to the band edges to s plane
  if (digital) {
    T <- 2       # sampling frequency of 2 Hz
    W <- 2 / T*tan(pi * W / T)
  }

  ## Generate splane poles and zeros for the chebyshev type 1 filter
  epsilon <- sqrt(10^(Rp/10) - 1)
  v0 <- asinh(1/epsilon)/n
  pole <- exp(1i*pi*c(seq(-(n-1), (n-1), by = 2))/(2*n))
  pole <- -sinh(v0)*Re(pole) + 1i*cosh(v0)*Im(pole)
  zero <- numeric(0)

  ## compensate for amplitude at s=0
  gain <- prod(-pole)
  ## if n is even, the ripple starts low, but if n is odd the ripple
  ## starts high. We must adjust the s=0 amplitude to compensate.
  if (n %% 2 == 0)
    gain <- gain/10^(Rp/20)

  ZPG <- Zpg(zero = zero, pole = pole, gain = gain)

  ## s-plane frequency transform
  ZPG <- sftrans(ZPG, W = W, stop = stop)

  ## Use bilinear transform to convert poles to the z plane
  if (digital)
     ZPG <- bilinear(ZPG, T = T)

  as.Arma(ZPG)
}
