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
##
## Completed by: Laurent S. Mazet

## Compute chebyshev type I filter order and cutoff for the desired response
## characteristics. Rp is the allowable decibels of ripple in the pass 
## band. Rs is the minimum attenuation in the stop band.
##
## [n, Wc] = cheb1ord(Wp, Ws, Rp, Rs)
##     Low pass (Wp<Ws) or high pass (Wp>Ws) filter design.  Wp is the
##     pass band edge and Ws is the stop band edge.  Frequencies are
##     normalized to [0,1], corresponding to the range [0,Fs/2].
## 
## [n, Wc] = cheb1ord([Wp1, Wp2], [Ws1, Ws2], Rp, Rs)
##     Band pass (Ws1<Wp1<Wp2<Ws2) or band reject (Wp1<Ws1<Ws2<Wp2)
##     filter design. Wp gives the edges of the pass band, and Ws gives
##     the edges of the stop band.
##
## See also: cheby1

cheb1ord <- function(Wp, Ws, Rp, Rs)  {

  if (length(Wp) != length(Ws))
    stop("Wp and Ws must have the same length")
  if (length(Wp) != 1 && length(Wp) != 2)
    stop("Wp,Ws must have length 1 or 2")
  if (length(Wp) == 2 && (all(Wp>Ws) || all(Ws>Wp) || diff(Wp)<=0 || diff(Ws)<=0))
    stop("Wp[1]<Ws[1]<Ws[2]<Wp[2] or Ws[1]<Wp[1]<Wp[2]<Ws[2]")

  T <- 2

  ## returned frequency is the same as the input frequency
  Wc <- Wp

  ## warp the target frequencies according to the bilinear transform
  Ws <- (2/T) * tan(pi *  Ws / T)
  Wp <- (2/T) * tan(pi * Wp / T)

  if (Wp[1] < Ws[1]) {
    ## low pass
    if (length(Wp) == 1) {
      Wa <- Ws/Wp
      type <- "low"
    } else {
      ## band reject
      type <- "stop"
      stop("band reject is not implement yet")
    }
  } else {
   ## if high pass, reverse the sense of the test
    if (length(Wp) == 1) {
      type <- "high"
      Wa <- Wp/Ws
    } else {
      type <- "pass"
      ## band pass 
      Wa<-(Ws^2 - Wp[1]*Wp[2]) / (Ws*(Wp[1]-Wp[2]))
    } 
  } 
  Wa <- min(abs(Wa))
  
  ## compute minimum n which satisfies all band edge conditions
  stop_atten <- 10^(abs(Rs)/10)
  pass_atten <- 10^(abs(Rp)/10)
  n <- ceiling(acosh(sqrt((stop_atten-1) / (pass_atten-1))) / acosh(Wa))
  FilterOfOrder(n = n, Wc = Wc, type = type, Rp = Rp)
} 
