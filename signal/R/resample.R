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

## usage: y=resample(x,p,q,d)
##
## Change the sample rate of x by a factor of p/q.  Note that p and q do
## not need to be integers since this routine does not use a polyphase
## rate change algorithm, but instead uses bandlimited interpolation,
## wherein the continous time signal is estimated by summing the sinc
## functions of the nearest neighbouring points up to distance d.
##
## This is discussed in:
##     J. O. Smith and P. Gossett (1984). A flexible sampling-rate
##     conversion method. In ICASSP-84, Volume II, pp. 19.4.1-19.4.2. 
##     New York: IEEE Press.
## See the authors page at: http://www-ccrma.stanford.edu/~jos/resample/
##
## Note that the resampling is not yet very fast or very good, but it is
## very flexible.
##
## Example
##    ## Speech example
##    [x, fs] = auload(file_in_loadpath("sample.wav"))
##    sound(resample(x,16000,fs), 16000);  # resample at 16 kHz
##
##    ## Example from interp1
##    xf=0:0.05:10.95; yf = sin(2*pi*xf/5)
##    xp=0:10;         yp = sin(2*pi*xp/5)
##    r = resample(yp,xp(2),xf(2))
##    plot(xf,yf,';original;',xf,r,';resample;',xp,yp,'*;;')
##
## Note that resample computes all samples up to but not including time
## n+1. If you are increasing the sample rate, this means that it will
## generate samples beyond the end of the time range of the original
## signal. That is why xf must goes all the way to 10.95 in the example.
 
## TODO: Fix so that audible clicking goes away.
## TODO: Change to a faster algorithm.
## TODO: Test on a chirp signal.
   
resample <- function(x, p, q = 1, d = 5) {
  order = d
  ## if rate reduction, apply antialiasing filter first
  r = p / q
  if (r < 1) {
    b = fir1(2*order+1, r)
    x = fftfilt(b, x)
  }

  ## Determine the new sampling times, and their distance to the old
  ## ones.  Note that the new series should be the maximum that can
  ## be contained in the old series without going over the time
  ## allotted to the old series.  In short, you have to go a little
  ## beyond the last sample of the old series if your new sampling
  ## rate is higher.
  t= seq(1, length(x)+1-1/r, by = 1/r)   # the sampling points of the new series
  idx = trunc(t)                         # the nearest old point
  t = t - idx                            # distance to the nearest old point

  ## generate the new series by summing the sinc functions of the
  ## nearest neighbour points implicit in the continuous time
  ## expansion of the old series.  This new series is truncated
  ## to +/- order nearest neighbours.  For convenience, the original
  ## series is zero-padded before and after, implicitly setting the
  ## neighbours at the start of the signal to zero.
  x = c(numeric(order), x, numeric(order))
  y = numeric(length(idx))        # the new series
  for (i in (-order):order) {
    w = sinc(t - i) * (0.5 + 0.5*cos(pi * (t-i) / (order + 0.5)))  # hanning window
    y = y + x[idx+i+order] * w
  }
  y
} 

#!demo
### xf=seq(0,10.95,by=0.05); yf = sin(2*pi*xf/5)
### xp=0:10;      yp = sin(2*pi*xp/5)
### r = resample(yp,xp[2],xf[2])
### title("confirm that the resampled function matches the original")
### plot(xf,yf, type = "l", col = "blue")
### lines(xf, r[1:length(xf)], col = "red")
### points(xp,yp, pch = 19, col = "blue")
### legend("bottomleft", c("Original", "Resample", "Data"),
###        col = c("blue", "red", "blue"),
###        pch = c(NA, NA, 19),
###        lty = c(1, 1, NA), bty = "n")
