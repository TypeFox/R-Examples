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

## Compute butterworth filter order and cutoff for the desired response
## characteristics. Rp is the allowable decibels of ripple in the pass 
## band. Rs is the minimum attenuation in the stop band.
##
## [n, Wc] = buttord(Wp, Ws, Rp, Rs)
##     Low pass (Wp<Ws) or high pass (Wp>Ws) filter design.  Wp is the
##     pass band edge and Ws is the stop band edge.  Frequencies are
##     normalized to [0,1], corresponding to the range [0,Fs/2].
## 
## [n, Wc] = buttord([Wp1, Wp2], [Ws1, Ws2], Rp, Rs)
##     Band pass (Ws1<Wp1<Wp2<Ws2) or band reject (Wp1<Ws1<Ws2<Wp2)
##     filter design. Wp gives the edges of the pass band, and Ws gives
##     the edges of the stop band.
##
## Theory: |H(W)|^2 = 1/[1+(W/Wc)^(2N)] = 10^(-R/10)
## With some algebra, you can solve simultaneously for Wc and N given
## Ws,Rs and Wp,Rp.  For high pass filters, subtracting the band edges
## from Fs/2, performing the test, and swapping the resulting Wc back
## works beautifully.  For bandpass and bandstop filters this process
## significantly overdesigns.  Artificially dividing N by 2 in this case
## helps a lot, but it still overdesigns.
##
## See also: butter

buttord <- function(Wp, Ws, Rp, Rs)  { 
  if (length(Wp) != length(Ws))
    stop("Wp and Ws must have the same length")
  if (length(Wp) != 1 && length(Wp) != 2)
    stop("Wp, Ws must have length 1 or 2")
  if (length(Wp) == 2 && (all(Wp>Ws) || all(Ws>Wp) || diff(Wp)<=0 || diff(Ws)<=0))
    stop("Wp(1)<Ws(1)<Ws(2)<Wp(2) or Ws(1)<Wp(1)<Wp(2)<Ws(2)")

  T <- 2
  
  ## if high pass, reverse the sense of the test
  stop <- which(Wp > Ws)
  Wp[stop] <- 1 - Wp[stop]   # stop will be at most length 1, so no need to
  Ws[stop] <- 1 - Ws[stop]   # subtract from matrix(1, 1,length(stop))

  if (length(Wp) == 2) {
    warning("buttord seems to overdesign bandpass and bandreject filters")
    type <- if (any(stop)) "stop" else "pass"
  } else {
    type <- if (any(stop)) "high" else "low"
  }

## Not sure why this was needed, but it generated wrong results:
#  if (any(stop)) type <- ""
## Quoting Andy Babour:
# I'm writing to see if 'signal::buttord'  is producing correct results.
# Here's an example (using version 0.7-3):
#
#    b <- buttord(.003, .001, 0.5, 29)
#    try(plot(freqz(butter(b)))) # error
#
#
#I don't doubt the value of the resulting filter order, but it sets
#'type' incorrectly to an empty string; whereas, the documentation
#implies it should be the string "high".  Conversely, if I flip Ws and Wp
#the result correctly shows type="low".
#
#I can easily get around this with
#
#    b$type <- "high"
#    plot(freqz(butter(b))) # ok 


  ## warp the target frequencies according to the bilinear transform
  Ws <- (2/T) * tan(pi * Ws / T)
  Wp <- (2/T) * tan(pi * Wp / T)
  
  ## compute minimum n which satisfies all band edge conditions
  ## the factor 1/length(Wp) is an artificial correction for the
  ## band pass/stop case, which otherwise significantly overdesigns.
  qs <- log(10^(Rs/10) - 1)
  qp <- log(10^(Rp/10) - 1)
  n <- ceiling(max(0.5*(qs - qp) / log(Ws/Wp)) /length(Wp))

  ## compute -3dB cutoff given Wp, Rp and n
  Wc <- exp(log(Wp) - qp/2/n)

  ## unwarp the returned frequency
  Wc <- atan(T/2*Wc)*T/pi
  
  ## if high pass, reverse the sense of the test
  Wc[stop] <- 1 - Wc[stop]
  FilterOfOrder(n = n, Wc = Wc, type = type)
} 
