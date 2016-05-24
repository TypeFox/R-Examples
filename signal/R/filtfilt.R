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

## usage: y = filtfilt(b, a, x)
##
## Forward and reverse filter the signal. This corrects for phase
## distortion introduced by a one-pass filter, though it does square the
## magnitude response in the process. That's the theory at least.  In
## practice the phase correction is not perfect, and magnitude response
## is distorted, particularly in the stop band.
##
## In this version, I zero-pad the end of the signal to give the reverse
## filter time to ramp up to the level at the end of the signal.
## Unfortunately, the degree of padding required is dependent on the
## nature of the filter and not just its order, so this function needs
## some work yet.
##
## Example
##    [b, a]=butter(3, 0.1);                   % 10 Hz low-pass filter
##    t = 0:0.01:1.0;                         % 1 second sample
##    x=sin(2*pi*t*2.3)+0.25*randn(size(t));  % 2.3 Hz sinusoid+noise
##    y = filtfilt(b,a,x); z = filter(b,a,x); % apply filter
##    plot(t,x,';data;',t,y,';filtfilt;',t,z,';filter;')

## Changelog:
## 2000 02 pkienzle@kienzle.powernet.co.uk
##      - pad with zeros to load up the state vector on filter reverse.
##      - add example

## TODO: In Matlab filtfilt `reduces filter startup transients by carefully
## TODO:    choosing initial conditions, and by prepending onto the input
## TODO:    sequence a short, reflected piece of the input sequence'.
## TODO:    Once filtic is written, use that here.
## TODO: My version seems to have similar quality to matlab, but both are
## TODO:    pretty bad.  They do remove gross lag errors, though.
## TODO: Note that if x is really long, it might be worth doing
## TODO:   the zero padding as a separate call to filter so that the
## TODO:   vector never has to be copied. E.g.,
## TODO:      [y, state] = filter(b,a,x)
## TODO:      tail = filter(b,a,matrix(0, 1,max(length(b),length(a))),state)
## TODO:      [tail, state] = filter(b,a,flipXX(tail))
## TODO:      y = flipXX(filter(b,a,flipXX(y), state))
## TODO:   Don't know for what n this would be faster, if any, but the
## TODO:   memory saving might be nice.



filtfilt <- function(filt, ...) UseMethod("filtfilt")

filtfilt.default <- function(filt, a, x, ...)  { 
    y = filter(filt, a, c(x, numeric(2 * max(length(a), length(filt)))))
    y = rev(filter(filt, a, rev(y)))[seq_along(x)]
    y
} 

filtfilt.Arma <- function(filt, x, ...) # IIR
  filtfilt(filt$b, filt$a, x)

filtfilt.Ma <- function(filt, x, ...) # FIR
  filtfilt(as.Arma(filt), x)

filtfilt.Zpg <- function(filt, x, ...) # Zero-pole-gain ARMA representation
  filtfilt(as.Arma(filt), x)
