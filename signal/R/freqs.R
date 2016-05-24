# Copyright (C) 2003 Julius O. Smith III
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# Octave is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Octave; see the file COPYING.  If not, write to the Free
# Software Foundation, 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# Usage: H = freqs(B,A,W)
#
# Compute the s-plane frequency response of the IIR filter B(s)/A(s) as 
# H = polyval(B,j*W) / polyval(A,j*W).  If called with no output
# argument, a plot of magnitude and phase are displayed.
#
# Example:
#   B = [1 2]; A = [1 1]
#   w = linspace(0,4,128)
#   freqs(B,A,w)

# 2003-05-16 Julius Smith - created

freqs <- function(filt, ...) UseMethod("freqs")

freqs.Arma <- function(filt, ...) # IIR
  freqs(filt$b, filt$a, ...)

freqs.Ma <- function(filt, ...) # FIR
  freqs.default(filt, 1, ...)

freqs.default <- function(filt = 1, a = 1, W, ...)  {
  
  H = polyval(filt,1i*W) / polyval(a,1i*W)
  
  res = list(H = H, W = W)
  class(res) = "freqs"
  res
} 

print.freqs <- plot.freqs <- function(x, ...)
  freqs_plot(x$W, x$H, ...)

freqs_plot <- function(w, ...) UseMethod("freqs_plot")

freqs_plot.freqs <- function(w, ...) 
  freqs_plot(w$W, w$H, ...)

freqs_plot.default <- function(w, h, ...) {
  mag = 20*log10(abs(h))
  phase = unwrap(Arg(h))
  op = par(mfrow=c(2,1), mar=c(4,4,2,1))
  plot(w, mag, xlab = "", ylab = "Magnitude (dB)", type = "l", ...)
  title('Frequency response plot by freqs')
  plot(w, phase/(2*pi), xlab = "Frequency (rad/sec)", ylab = "Cycles", type = "l", ...)
  par(op)
}

#!demo
##B = c(1, 2)
##A = c(1, 1)
##w = seq(0,4,length=128)
##freqs(B,A,w)
