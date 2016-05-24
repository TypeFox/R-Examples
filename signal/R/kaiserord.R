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

## usage: [n, Wn, beta, ftype] = kaiserord(f, m, dev [, fs])
##
## Returns the parameters needed for fir1 to produce a filter of the
## desired specification from a kaiser window:
##       n: order of the filter (length of filter minus 1)
##       Wn: band edges for use in fir1
##       beta: parameter for kaiser window of length n+1
##       ftype: choose between pass and stop bands
##       b = fir1(n,Wn,kaiser(n+1,beta),ftype,'noscale')
##
## f: frequency bands, given as pairs, with the first half of the
##    first pair assumed to start at 0 and the last half of the last
##    pair assumed to end at 1.  It is important to separate the
##    band edges, since narrow transition regions require large order
##    filters.
## m: magnitude within each band.  Should be non-zero for pass band
##    and zero for stop band.  All passbands must have the same
##    magnitude, or you will get the error that pass and stop bands
##    must be strictly alternating.
## dev: deviation within each band.  Since all bands in the resulting
##    filter have the same deviation, only the minimum deviation is
##    used.  In this version, a single scalar will work just as well.
## fs: sampling rate.  Used to convert the frequency specification into
##    the [0, 1], where 1 corresponds to the Nyquist frequency, fs/2.
##
## The Kaiser window parameters n and beta are computed from the
## relation between ripple (A=-20*log10(dev)) and transition width 
## (dw in radians) discovered empirically by Kaiser:
##
##           / 0.1102(A-8.7)                        A > 50
##    beta = | 0.5842(A-21)^0.4 + 0.07886(A-21)     21 <= A <= 50
##           \ 0.0                                  A < 21
##
##    n = (A-8)/(2.285 dw)
##
## Example
##    [n, w, beta, ftype] = kaiserord([1000,1200], [1,0], [0.05,0.05], 11025)
##    freqz(fir1(n,w,kaiser(n+1,beta),ftype,'noscale'),1,[],11025)

## TODO: order is underestimated for the final test case: 2 stop bands.
## TODO:     octave> ftest("kaiserord") # shows test cases

kaiserord <- function(f, m, dev, Fs = 2)  { 

  ## parameter checking
  if (length(f) != 2*length(m)-2)
    stop("kaiserord must have one magnitude for each frequency band")
  if (length(m) > 2 && any(m[1:(length(m)-2)] != m[3:length(m)]))
    stop("kaiserord pass and stop bands must be strictly alternating")
  if (length(dev) != length(m) && length(dev) != 1)
    stop("kaiserord must have one deviation for each frequency band")
  dev = min(dev)
  if (dev <= 0)
    stop("kaiserord must have dev > 0")

  ## use midpoints of the transition region for band edges
  w = (f[seq(1, length(f), by = 2)] + f[seq(2, length(f), by = 2)])/Fs

  ## determine ftype
  if (length(w) == 1)
    if (m[1] > m[2])  ftype='low' else ftype='high'
  else if (length(w) == 2)
    if (m[1] > m[2]) ftype='stop' else ftype='pass'
  else
    if (m[1] > m[2]) ftype='DC-1' else ftype='DC-0'

  ## compute beta from dev
  A = -20*log10(dev)
  if (A > 50)
    beta = 0.1102*(A-8.7)
  else if (A >= 21)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21)
  else
    beta = 0.0

  ## compute n from beta and dev
  dw = 2*pi*min(f[seq(2, length(f), by = 2)] - f[seq(1, length(f), by = 2)])/Fs
  n = max(1, ceiling((A-8) / (2.285*dw)))

  ## if last band is high, make sure the order of the filter is even.
  if ((m[1] > m[2]) == (length(w) %% 2 == 0) &&
      n %% 2 == 1)
    n = n+1

  FilterOfOrder(n = n, Wc = w, type = ftype, beta = beta)
} 

#!demo
###Fs = 11025
###op = par(mfrow=c(2,2), mar=c(3,3,1,1))
###for ( i in 1 : 4) {
###  if (i==1) {
###    bands=c(1200, 1500); mag=c(1, 0); dev=c(0.1, 0.1)
###  } else if (i==2) {
###    bands=c(1000, 1500); mag=c(0, 1); dev=c(0.1, 0.1)
###  } else if (i==3) {
###    bands=c(1000, 1200, 3000, 3500); mag=c(0, 1, 0); dev=0.1
###  } else if (i==4) {
###    bands=100*c(10, 13, 15, 20, 30, 33, 35, 40)
###    mag=c(1, 0, 1, 0, 1); dev=0.05
###  }
###  kaisprm = kaiserord(bands, mag, dev, Fs)
###  with(kaisprm, {
###    d <<- max(1,trunc(n/10))
###    if (mag[length(mag)]==1 && d %% 2 ==1)
###      d<<-d+1
###    f1 <<- freqz(fir1(n,w,ftype,kaiser(n+1,beta),'noscale'), Fs = Fs, plot = FALSE)
###    f2 <<- freqz(fir1(n-d,w,ftype,kaiser(n-d+1,beta),'noscale'), Fs = Fs, plot = FALSE)
###  })                                                               
###  plot(f1$f,abs(f1$h), col = "blue", type="l", xlab = "", ylab = "")
###  lines(f2$f,abs(f2$h), col = "red")
###  legend("right", paste("order", c(kaisprm$n-d, kaisprm$n)), col=c("red", "blue"), lty=1, bty="n")
###  b = c(0, bands, Fs/2)
###  for ( i in seq(2,length(b),by=2)) {
###    hi=mag[i/2]+dev[1]; lo=max(mag[i/2]-dev[1],0)
###    lines(c(b[i-1], b[i], b[i], b[i-1], b[i-1]),c(hi, hi, lo, lo, hi))
###  } 
###}
###par(op)
#! %--------------------------------------------------------------
#! % A filter meets the specifications if its frequency response
#! % passes through the ends of the criteria boxes, and fails if
#! % it passes through the top or the bottom.  The criteria are
#! % met precisely if the frequency response only passes through
#! % the corners of the boxes.  The blue line is the filter order
#! % returned by kaiserord, and the red line is some lower filter
#! % order.  Confirm that the blue filter meets the criteria and
#! % the red line fails.

#!# XXX FIXME XXX end demo to show detail at criteria box corners

