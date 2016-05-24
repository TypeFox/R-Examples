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

## usage: [x, t] = impz(b [, a, n, fs])
##
## Generate impulse-response characteristics of the filter. The filter
## coefficients correspond to the the z-plane rational function with
## numerator b and denominator a.  If a is not specified, it defaults to
## 1. If n is not specified, or specified as [], it will be chosen such
## that the signal has a chance to die down to -120dB, or to not explode
## beyond 120dB, or to show five periods if there is no significant
## damping. If no return arguments are requested, plot the results.
##
## See also: freqz, zplane

## 1999 pkienzle@kienzle.powernet.co.uk
##      - if nargout=0, produce plot and don't set return values

## TODO: Call equivalent function from control toolbox since it is
## TODO:    probably more sophisticated than this one, and since it
## TODO:    is silly to maintain two different versions of essentially
## TODO:    the same thing.

impz <- function(filt, ...) UseMethod("impz")

print.impz <- plot.impz <- function(x, xlab = "Time, msec", ylab = "", type = "l",
    main = "Impulse response", ...) {
  plot(x$t, x$x, xlab = xlab, ylab = ylab, type = type, main = main, ...)
}

impz.Arma <- function(filt, ...) # IIR
  impz(filt$b, filt$a, ...)

impz.Ma <- function(filt, ...) # FIR
  impz(filt$b, 1, ...)

impz.default <- function(filt, a = 1, n = NULL, Fs = 1, ...)  {
  b = filt
  if (length(n) == 0 && length(a) > 1) {
    precision = 1e-6
    r = roots(a)
    maxpole = max(abs(r))
    if (maxpole > 1+precision)     # unstable -- cutoff at 120 dB
      n = floor(6/log10(maxpole))
    else if (maxpole < 1-precision) # stable -- cutoff at -120 dB
      n = floor(-6/log10(maxpole))
    else {                       # periodic -- cutoff after 5 cycles
      n = 30
      ## find longest period less than infinity
      ## cutoff after 5 cycles (w=10*pi)
      rperiodic = r[abs(r) >= 1 - precision & abs(Arg(r)) > 0]
      if (!is.null(rperiodic)) {
        n_periodic = ceiling(10*pi / min(abs(Arg(rperiodic))))
        if (n_periodic > n)
          n = n_periodic
      } 
      
      ## find most damped pole
      ## cutoff at -60 dB
      rdamped = r[abs(r) < 1 - precision]
      if (!is.null(rdamped))
      n_damped = floor(-3/log10(max(abs(rdamped))))
      if (n_damped > n)
        n = n_damped
    } 
    n = n + length(b)
  } else if (is.null(n)) {
    n = length(b)
  }
  if (length(a) == 1)
    x = fftfilt(b/a, c(1, numeric(n-1)))
  else
    x = filter(b, a, c(1, numeric(n-1)))
  t = (0:(n-1))/Fs
  
  res = list(x = x, t = t)
  class(res) = "impz"
  res
} 
