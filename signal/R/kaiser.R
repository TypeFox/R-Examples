## Copyright (C) 1995, 1996, 1997  Kurt Hornik
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details. 
## 
## You should have received a copy of the GNU General Public License
## along with this file.  If not, write to the Free Software Foundation,
## 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

## usage:  kaiser (n, beta)
##
## Returns the filter coefficients of the n-point Kaiser window with
## parameter beta.
##
## For the definition of the Kaiser window, see A. V. Oppenheim &
## R. W. Schafer, "Discrete-Time Signal Processing".
##
## The continuous version of width n centered about x=0 is:
##
##         besseli(0, beta * sqrt(1-(2*x/n).^2))
## k(x) =  -------------------------------------,  n/2 <= x <= n/2
##                besseli(0, beta)
##
## See also: kaiserord
  
## Author:  Kurt Hornik <Kurt.Hornik@ci.tuwien.ac.at>
## Description:  Coefficients of the Kaiser window

## 2000-02 Paul Kienzle (pkienzle@kienzle.powernet.co.uk)
##    use besseli rather than jybess
##    note, although Oppenheim & Schafer, 2nd edition has a formula
##    which looks completely different than the one herein, it gives
##    identical results
  
kaiser <- function(n, beta)  { 
  
  if ( !(length(n) == 1 && (n == round(n)) && (n > 0))) 
    stop("kaiser:  n has to be a positive integer")

  if ( !(length(beta) == 1 && (beta == as.double(beta))))
    stop("kaiser:  beta has to be a real scalar")
  
  if (n == 1)
    w = 1
  else {
    m = n - 1
    k = 0:m
    k = 2 * beta / m * sqrt (k * (m - k))
    w = besselI(k, 0) / besselI(beta, 0)
  }
  w
}
