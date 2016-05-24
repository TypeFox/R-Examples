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
## @deftypefn {Function File} {[@var{h}, @var{w}] =} freqz (@var{b}, @var{a}, @var{n}, "whole")
## Return the complex frequency response @var{h} of the rational IIR filter
## whose numerator and denominator coefficients are @var{b} and @var{a},
## respectively.  The response is evaluated at @var{n} angular frequencies
## between 0 and
## @ifinfo
##  2*pi.
## @} # ifinfo
## @iftex
## @tex
##  $2\pi$.
## @} # tex
## @} # iftex
##
## @noindent
## The output value @var{w} is a vector of the frequencies.
##
## If the fourth argument is omitted, the response is evaluated at
## frequencies between 0 and
## @ifinfo
##  pi.
## @} # ifinfo
## @iftex
## @tex
##  $\pi$.
## @} # tex
## @} # iftex
##
## If @var{n} is omitted, a value of 512 is assumed.
##
## If @var{a} is omitted, the denominator is assumed to be 1 (this
## corresponds to a simple FIR filter).
##
## For fastest computation, @var{n} should factor into a small number of
## small primes.
##
## @deftypefnx {Function File} {@var{h} =} freqz (@var{b}, @var{a}, @var{w})
## Evaluate the response at the specific frequencies in the vector @var{w}.
## The values for @var{w} are measured in radians.
##
## @deftypefnx {Function File} {[@dots{}] =} freqz (@dots{}, @var{Fs})
## Return frequencies in Hz instead of radians assuming a sampling rate
## @var{Fs}.  If you are evaluating the response at specific frequencies 
## @var{w}, those frequencies should be requested in Hz rather than radians.
##
## @deftypefnx {Function File} {} freqz (@dots{})
## Plot the pass band, stop band and phase response of @var{h} rather
## than returning them.
## @} # deftypefn

## Author: jwe ???

freqz <- function(filt, ...) UseMethod("freqz")

freqz.Arma <- function(filt, ...) # IIR
  freqz(filt$b, filt$a, ...)

freqz.Ma <- function(filt, ...) # FIR
  freqz.default(filt, 1, ...)

print.freqz <- plot.freqz <- function(x, ...)
  freqz_plot(x$f, x$h)

freqz.default <- function(filt = 1, a = 1, n = 512, region = NULL, Fs = 2*pi, ...)  {
  b <- filt
  if (is.null(region))
    region <- if (is.double(b) && is.double(a)) "half" else "whole"
  if (length(n) > 1) { ## Explicit frequency vector given
    f <- n
    w <- 2*pi*f/Fs
    hb <- polyval(rev(b), exp(-1i*w))
    ha <- polyval(rev(a), exp(-1i*w))
  } else if (region == "whole") {
    f <- Fs * (0:(n-1)) / n
    ## polyval(fliplr(P),exp(-jw)) is O(p n) and fft(x) is O(n log(n)), where p is the 
    ## order of the the polynomial P.  For small p it would be faster to use polyval  
    ## but in practice the overhead for polyval is much higher and the little bit of
    ## time saved isn't worth the extra code.
    hb <- fft(postpad(b, n))
    ha <- fft(postpad(a, n))
  } else { # region == "half"
    f <- Fs/2 * (0:(n-1)) / n
    hb <- fft(postpad(b, 2*n))[1:n]
    ha <- fft(postpad(a, 2*n))[1:n]
  }

  h <- hb / ha

  res <- list(h = h, f = f)
  class(res) <- "freqz"
  res
} 

#!test # correct values and fft-polyval consistency
#! # butterworth filter, order 2, cutoff pi/2 radians
#! b = [0.292893218813452  0.585786437626905  0.292893218813452]
#! a = [1  0  0.171572875253810]
#! [h,w] = freqz(b,a,32)
#! assert(h(1),1,10*eps)
#! assert(abs(h(17)).^2,0.5,10*eps)
#! assert(h,freqz(b,a,w),10*eps); # fft should be consistent with polyval

#!test # whole-half consistency
#! b = c(1,1,1)/3; # 3-sample average
#! res = freqz(b,1,32,'whole',plot=F)
#! assert(h(2:16),conj(h(32:-1:18)),20*eps)
#! [h2,w2] = freqz(b,1,16,'half')
#! assert(h(1:16),h2,20*eps)
#! assert(w(1:16),w2,20*eps)

#!test # Sampling frequency properly interpreted
#! b = [1 1 1]/3
#! [h,f] = freqz(b,1,16,320)
#! assert(f,[0:15]*10,10*eps)
#! [h2,f2] = freqz(b,1,[0:15]*10,320)
#! assert(f2,[0:15]*10,10*eps)
#! assert(h,h2,20*eps)
#! [h3,f3] = freqz(b,1,32,'whole',320)
#! assert(f3,[0:31]*10,10*eps)

## Based on code originally by John W. Eaton
## Converted to R by Tom Short, tshort@eprisolutions.com
## Changes Copyright 2006 EPRI Solutions, Inc. also under the GNU GPL
## version 2 or later.

## Copyright (C) 2002 John W. Eaton
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
## @deftypefn {Function File} freqz_plot (@var{w}, @var{h})
## Plot the pass band, stop band and phase response of @var{h}.
## @} # deftypefn

## Author: Paul Kienzle <pkienzle@users.sf.net>

freqz_plot <- function(w, ...) UseMethod("freqz_plot")

freqz_plot.freqz <- function(w, ...) 
  freqz(w$h, w$f, ...)

freqz_plot.default <- function(w, h, ...) {
  mag = 20 * log10(abs(h))
  phase = unwrap(Arg(h))
  maxmag = max(mag)
  op = par(mfrow=c(3,1), mar = c(4,4,1.5,1))
  on.exit(par(op))
  plot(w, mag, type = "l", xlab = "", ylab = "", ylim = c(maxmag-3, maxmag), ...)
  title("Pass band (dB)")
  plot(w, mag, type = "l", xlab = "", ylab = "", ...)
  title("Stop band (dB)")
  plot(w, phase*360/(2*pi), type = "l", ..., xlab = "Frequency", ylab = "")
  title("Phase (degrees)")
} 
