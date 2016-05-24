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

zplane <- function(filt, ...) UseMethod("zplane")

zplane.Arma <- function(filt, ...) # IIR
  zplane(filt$b, filt$a, ...)

zplane.Ma <- function(filt, ...) # FIR
  zplane(filt$b, 1, ...)

zplane.Zpg <- function(filt, ...) {
  x = filt
  r = exp(2i*pi*(0:100)/100)
  xlim = range(c(-1.1, 1.1, Re(x$pole), Re(x$zero)))
  ylim = range(c(-1.1, 1.1, Im(x$pole), Im(x$zero)))
  plot(Re(r), Im(r), col = "red", xlab = "", ylab = "", xlim = xlim, ylim = ylim, type = "l", asp = 1, ...)
  points(Re(x$pole), Im(x$pole), pch = 4)
  points(Re(x$zero), Im(x$zero))
}

zplane.default <- function(filt, a, ...) {
  zplane(Zpg(roots(filt), roots(a), 1), ...)
}

#!demo
#! ## construct target system:
#! ##   symmetric zero-pole pairs at r*exp(iw),r*exp(-iw)
#! ##   zero-pole singletons at s
#! pw=[0.2, 0.4, 0.45, 0.95];   #pw = [0.4]
#! pr=[0.98, 0.98, 0.98, 0.96]; #pr = [0.85]
#! ps=[]
#! zw=[0.3];  # zw=[]
#! zr=[0.95]; # zr=[]
#! zs=[]
#! 
#! ## system function for target system
#! p=[[pr, pr].*exp(1i*pi*[pw, -pw]), ps]'
#! z=[[zr, zr].*exp(1i*pi*[zw, -zw]), zs]'
#! M = length(z); N = length(p)
#! sys_a = [ matrix(0, 1, M-N), real(poly(p)) ]
#! sys_b = [ matrix(0, 1, N-M), real(poly(z)) ]

#! save_replot = automatic_replot
#! automatic_replot = 0
#! disp("The first two graphs should be identical, with poles at (r,w)=")
#! disp(sprintf(" (%.2f,%.2f)", [pr ; pw]))
#! disp("and zeros at (r,w)=")
#! disp(sprintf(" (%.2f,%.2f)", [zr ; zw]))
#! disp("with reflection across the horizontal plane")
#! subplot(231); title("transfer function form"); zplane(sys_b, sys_a)
#! subplot(232); title("pole-zero form"); zplane(z,p)

#! subplot(233); title("empty p"); zplane(z)
#! subplot(234); title("empty a"); zplane(sys_b)
#! disp("The matrix plot has 2 sets of points, one inside the other")
#! subplot(235); title("matrix"); zplane([z, 0.7*z], [p, 0.7*p])
#! oneplot()
#! automatic_replot = save_replot
