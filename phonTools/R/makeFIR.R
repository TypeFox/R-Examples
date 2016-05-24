# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


makeFIR = function (frequency, dB, order = 200, signal = NULL, window = 'hann', verify = FALSE, interpolation = 'linear'){
  frequency = sort(frequency)  
  if (order %% 2) order = order + 1
  
  if ((frequency[1]) != 0) stop ('First frequency point must be 0.')
  if (length(dB) != length (frequency)) stop ('Frequenecy and dB specifications of unequal length.')
  if (order < 10) stop ('Filter order must be at least 10 taps.')
  
  n = length (frequency)
  if (n < 2) stop ('At least two points are required to specify the filter.')
  
  dB = dB - max (dB)
  ylim = range (dB); ylim[1] = ylim[1]-20; ylim[2] = ylim[2]+5;
  
  spots = seq (0, tail(frequency,1), length.out = order/2 + 1)
  for (i in 1:length (frequency)) frequency[i] = spots[order (abs(spots - frequency[i]))[1]]
  
  incr = spots[2] - spots[1]  
  dBs = interpolate (dB, frequency, steps = 20, increment = incr, show = F, type = interpolation)[,2]
  l = length (dBs)
  dBs = c(dBs, rev(dBs[-c(1,l)]))
  
  dBs = 10^(dBs/20)
  response = Re(fft (dBs, inverse = TRUE))
  response = c(response[l:1], rev(response[l:1][-l])) 
  response = (response / max(response)) * windowfunc(response, type = window)
  
  if (!is.null(signal)) output = FIRfilter (signal, impulse = response) 
  
  oldpar = par()
  if (verify){
    oldpar = par()
	if(is.null(signal)){
		par (mfrow = c(2,1), mar = c(4.5,4.5,3,1))
		plot (response, main = 'Filter Impulse Response', xlab = 'Tap', 
			  ylab = 'Amplitude', xaxs='i', type = 'b')
		spectralslice (response, fs = max(frequency)*2, padding = 15000, main = 'Filter Frequency Response', 
					   ylim = ylim, window = window)
		points (frequency, dB, cex = 1.3, col = 4, pch = 16)
	}  
    if (verify & !is.null(signal)){
		par (mfrow = c(2,2), mar = c(4.5,4.5,3,1))
		plot (response, main = 'Filter Impulse Response', xlab = 'Tap', 
			  ylab = 'Amplitude', xaxs='i', type = 'b')
		spectralslice (response, fs = max(frequency)*2, padding = 15000, main = 'Filter Frequency Response', 
					   ylim = ylim, window = window)
		points (frequency, dB, cex = 1.3, col = 4, pch = 16)
		spectralslice (signal, fs = max(frequency)*2, padding = 5000, main = 'Pre-Filtering', 
					   ylim = ylim, window = window)
		spectralslice (output, fs = max(frequency)*2, padding = 5000, main = 'Post-Filtering', 
					   ylim = ylim, window = window)
		spectralslice (response, fs = max(frequency)*2, padding = 15000, main = 'Filter Frequency Response', 
					   ylim = ylim, window = window, add = TRUE, lty = 'dotted', col = 4)
		par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1));
	  }
	suppressWarnings (par (oldpar)) 
  } 
  if (is.null(signal)) return (response)
  if (!is.null(signal)) return (signal)
}

