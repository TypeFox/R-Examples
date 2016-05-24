## Copyright (C) 2002 Andras Carezia
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
## USA

## Usage:  chebwin (n, at)
##
## Returns the filter coefficients of the n-point Dolph-Chebyshev window
## with at dB of attenuation in the stop-band of the corresponding
## Fourier transform.
##
## For the definition of the Chebyshev window, see
##
## * Peter Lynch, "The Dolph-Chebyshev Window: A Simple Optimal Filter",
##   Monthly Weather Review, Vol. 125, pp. 655-660, April 1997.
##   (http://www.maths.tcd.ie/~plynch/Publications/Dolph.pdf)
##
## * C. Dolph, "A current distribution for broadside arrays which
##   optimizes the relationship between beam width and side-lobe level",
##   Proc. IEEE, 34, pp. 335-348.
##
## The window is described in frequency domain by the expression:
##
##          Cheb(n-1, beta * cos(pi * k/n))
##   W(k) = -------------------------------
##                 Cheb(n-1, beta)
##
## with
##
##   beta = cosh(1/(n-1) * acosh(10^(at/20))
##
## and Cheb(m,x) denoting the m-th order Chebyshev polynomial calculated
## at the point x.
##
## Note that the denominator in W(k) above is not computed, and after
## the inverse Fourier transform the window is scaled by making its
## maximum value unitary.
##
## See also: kaiser

## $Id: chebwin.m,v 1.4 2005/12/29 03:54:39 pkienzle Exp $
##
## Author:  Andras Carezia <acarezia@uol.com.br>
## Description:  Coefficients of the Dolph-Chebyshev window

chebwin  <- function(n, at)  {

  if (!(length(n) == 1 && (n == round(n)) && (n > 0)))
    stop("n has to be a positive integer")
  if (!(length(at) == 1 && (at == Re(at))))
    stop("at has to be a real scalar")
  
  if (n == 1)
    w <- 1
  else {
    # beta calculation
    gamma <- 10^(-at/20)
    beta <- cosh(1/(n-1) * acosh(1/gamma))
    # freq. scale
    k <- 0:(n-1)
    x <- beta*cos(pi*k/n)
    # Chebyshev window (freq. domain)
    p <- cheb(n-1, x)
    # inverse Fourier transform
    if (n %% 2) {
      w <- Re(fft(p))
      M <- (n+1)/2
      w <- w[1:M] / w[1]
      w <- c(w[M:2], w)
    } else {
      # half-sample delay (even order)
      p <- p * exp(1i*pi/n * (0:(n-1)))
      w <- Re(fft(p))
      M <- n / 2 + 1
      w <- w/w[2]
      w <- c(w[M:2], w[2:M])
    }
  }
  w
}

## Copyright (C) 2002 Andras Carezia
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
## USA

## Usage:  cheb (n, x)
##
## Returns the value of the nth-order Chebyshev polynomial calculated at
## the point x. The Chebyshev polynomials are defined by the equations:
##
##           / cos(n acos(x),    |x| <= 1
##   Tn(x) = |
##           \ cosh(n acosh(x),  |x| > 1
##
## If x is a vector, the output is a vector of the same size, where each
## element is calculated as y(i) = Tn(x(i)).

## Author:  Andras Carezia <acarezia@uol.com.br>
## Description:  Value of the Chebyshev polynomials

cheb  <- function(n, x)  {
  if (!(is.numeric(n) && (n == round(n)) && (n >= 0)))
    stop("n has to be a positive integer")

  T <- numeric(length(x))
  
  ind <- x <= 1
  if (any(ind))
    T[ind] <- cos(n * acos(as.complex(x[ind])))

  ind <- x > 1
  myacosh <- function(x) log(x + sqrt(x^2 - 1)) # workaround for a win32 bug in acosh
  if (any(ind))
    T[ind] <- cosh(n * myacosh(as.complex(x[ind])))

  Re(T)
} 
