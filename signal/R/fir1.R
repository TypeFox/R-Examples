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

## usage: b = fir1(n, w [, type] [, window] [, noscale])
##
## Produce an order n FIR filter with the given frequency cutoff,
## returning the n+1 filter coefficients in b.  
##
## n: order of the filter (1 less than the length of the filter)
## w: band edges
##    strictly increasing vector in range [0, 1]
##    singleton for highpass or lowpass, vector pair for bandpass or
##    bandstop, or vector for alternating pass/stop filter.
## type: choose between pass and stop bands
##    'high' for highpass filter, cutoff at w
##    'stop' for bandstop filter, edges at w = [lo, hi]
##    'DC-0' for bandstop as first band of multiband filter
##    'DC-1' for bandpass as first band of multiband filter
## window: smoothing window
##    defaults to hamming(n+1) row vector
##    returned filter is the same shape as the smoothing window
## noscale: choose whether to normalize or not
##    'scale': set the magnitude of the center of the first passband to 1
##    'noscale': don't normalize
##
## To apply the filter, use the return vector b:
##       y=filter(b,1,x)
##
## Examples:
##   freqz(fir1(40,0.3))
##   freqz(fir1(15,[0.2, 0.5], 'stop'));  # note the zero-crossing at 0.1
##   freqz(fir1(15,[0.2, 0.5], 'stop', 'noscale'))

## TODO: Consider using exact expression (in terms of sinc) for the
## TODO:    impulse response rather than relying on fir2.
## TODO: Find reference to the requirement that order be even for
## TODO:    filters that end high.  Figure out what to do with the
## TODO:    window in these cases

fir1 <- function(n, w, type = c("low", "high", "stop", "pass", "DC-0", "DC-1"), 
    window = hamming(n+1), scale = TRUE)  { 

  type <- match.arg(type)
  if (!is.logical(scale)) {
    scale <- match.arg(scale, c("scale", "noscale"))
    scale <- scale == "scale"
  }
  if(is.function(window))
    window <- window(n+1)
  else if(is.character(window))
    window <- do.call(window, list(n+1))

  ## Assign default window, filter type and scale.
  ## If single band edge, the first band defaults to a pass band to 
  ## create a lowpass filter.  If multiple band edges, the first band 
  ## defaults to a stop band so that the two band case defaults to a 
  ## band pass filter.  Ick.
  ftype <- tolower(type) %in% c('low','stop','dc-1')
  
  ## build response function according to fir2 requirements
  bands <- length(w) + 1
  f <- numeric(2*bands)
  f[2*bands] = 1
  f[seq(2, 2*bands-1, by = 2)] <- w
  f[seq(3, 2*bands-1, by = 2)] <- w
  m <- numeric(2*bands)
  m[seq(1, 2*bands, by = 2)] <- (1:bands - (1-ftype)) %% 2
  m[seq(2, 2*bands, by = 2)] <- m[seq(1, 2*bands, by = 2)]

  ## Increment the order if the final band is a pass band.  Something
  ## about having a nyquist frequency of zero causing problems.
  if (n %% 2 == 1 && m[2*bands] == 1) { 
    warning("n must be even for highpass and bandstop filters. Incrementing.")
    n <- n + 1
    if (is.vector(window) && is.double(window)) {
      ## End the window using interpolation
      M <- length(window)
      if (M == 1)
        window <- c(window, window)
      else 
        window <- interp1(seq(0,1,length=M), window, seq(0,1,length=M+1), 
                          if (M < 4) 'linear' else 'spline')
    }
  }

  ## compute the filter
  b <- fir2(n, f, m, 512, 2, window)

  ## normalize filter magnitude
  if (scale) {
    ## find the middle of the first band edge
    if (m[1] == 1)
      w_o <- (f[2] - f[1])/2
    else
      w_o <- f[3] + (f[4] - f[3])/2

    ## compute |h(w_o)|^-1
    renorm <- 1/abs(polyval(b, exp(-1i*pi*w_o)))

    ## normalize the filter
    b <- renorm*b
  }
  Ma(b)
} 

#!demo
#! freqz(fir1(40,0.3))
#!demo
#! freqz(fir1(15,[0.2, 0.5], 'stop'));  # note the zero-crossing at 0.1
#!demo
#! freqz(fir1(15,[0.2, 0.5], 'stop', 'noscale'))

#!assert(fir1(2, .5, 'low', @hanning, 'scale'), [0 1 0]')
#!assert(fir1(2, .5, 'low', "hanning", 'scale'), [0 1 0]')
#!assert(fir1(2, .5, 'low', hamming(3), 'scale'), [0 1 0]')

#!assert(fir1(10,.5,'noscale'), fir1(10,.5,'low','hamming','noscale'))
#!assert(fir1(10,.5,'high'), fir1(10,.5,'high','hamming','scale'))
#!assert(fir1(10,.5,'boxcar'), fir1(10,.5,'low','boxcar','scale'))
#!assert(fir1(10,.5,'hanning','scale'), fir1(10,.5,'scale','hanning','low'))
#!assert(fir1(10,.5,'haNNing','NOscale'), fir1(10,.5,'noscale','Hanning','LOW'))
#!assert(fir1(10,.5,'boxcar',[]), fir1(10,.5,'boxcar'))
