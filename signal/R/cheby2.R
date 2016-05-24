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

## Generate an Chebyshev type II filter with Rs dB of pass band ripple.
## 
## [b, a] = cheby2(n, Rs, Wc)
##    low pass filter with cutoff pi*Wc radians
##
## [b, a] = cheby2(n, Rs, Wc, 'high')
##    high pass filter with cutoff pi*Wc radians
##
## [b, a] = cheby2(n, Rs, [Wl, Wh])
##    band pass filter with edges pi*Wl and pi*Wh radians
##
## [b, a] = cheby2(n, Rs, [Wl, Wh], 'stop')
##    band reject filter with edges pi*Wl and pi*Wh radians
##
## [z, p, g] = cheby2(...)
##    return filter as zero-pole-gain rather than coefficients of the
##    numerator and denominator polynomials.
##
## [...] = cheby2(...,'s')
##     return a Laplace space filter, W can be larger than 1.
## 
## [a,b,c,d] = cheby2(...)
##  return  state-space matrices 
## 
## References: 
##
## Parks & Burrus (1987). Digital Filter Design. New York:
## John Wiley & Sons, Inc.

## Author: Paul Kienzle <pkienzle@users.sf.net>
## Modified: Doug Stewart Feb. 2003

cheby2 <- function(n, ...) UseMethod("cheby2")

cheby2.FilterOfOrder <- function(n, ...)
  cheby1(n$n, n$Rp, n$Wc, n$type, ...)

cheby2.default <- function(n, Rp, W, type = c("low", "high", "stop", "pass"), plane = c("z", "s"), ...) {

  type <- match.arg(type)
  plane <- match.arg(plane)

  ## interpret the input parameters
  if (!(length(n)==1 && n == round(n) && n > 0))
    stop("cheby2: filter order n must be a positive integer")

  stop <- type == "stop" || type == "high"
  digital <- plane == "z"

  if (length(W) != 1 && length(W) != 2)
    stop("cheby2: frequency must be given as w0 or c(w0, w1)")

  if (digital && !all(W >= 0 & W <= 1))
    stop("cheby2: critical frequencies must be in (0 1)")
  else if (!digital && !all(W >= 0))
    stop("cheby2: critical frequencies must be in (0 inf)")

  if (Rp < 0)
    stop("cheby2: passband ripple must be positive decibels")

  ## Prewarp to the band edges to s plane
  if (digital) {
    T <- 2       # sampling frequency of 2 Hz
    W <- 2 / T*tan(pi * W / T)
  }

  ## Generate splane poles and zeros for the chebyshev type 2 filter
  ## From: Stearns, SD; David, RA; (1988). Signal Processing Algorithms. 
  ##       New Jersey: Prentice-Hall.
  C <- 1             # default cutoff frequency
  lambda <- 10^(Rp/20)
  phi <- log(lambda + sqrt(lambda^2-1))/n
  theta <- pi*((1:n) - 0.5)/n
  alpha <- -sinh(phi)*sin(theta)
  beta <- cosh(phi)*cos(theta)
  if (n %% 2)
    ## drop theta==pi/2 since it results in a zero at infinity
    zero <- 1i * C / cos(theta[c(1:((n-1)/2), ((n+3)/2):n)])
  else
   zero <- 1i * C / cos(theta)
  pole <- C / (alpha^2 + beta^2) * (alpha - 1i*beta)

  ## Compensate for amplitude at s=0
  ## Because of the vagaries of floating point computations, the
  ## prod(pole)/prod(zero) sometimes comes out as negative and
  ## with a small imaginary component even though analytically
  ## the gain will always be positive, hence the abs(Re(...))
  gain <- abs(Re(prod(pole) / prod(zero)))

  ZPG <- Zpg(zero = zero, pole = pole, gain = gain)

  ## s-plane frequency transform
  ZPG <- sftrans(ZPG, W = W, stop = stop)

  ## Use bilinear transform to convert poles to the z plane
  if (digital)
     ZPG <- bilinear(ZPG, T = T)

  as.Arma(ZPG)
}
