# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

sinusoid = function (freqs, amps = rep(1, length(freqs)), dur = 50, phases = rep(0, length(freqs)), fs = 10000, sum = FALSE, show = FALSE, colors = NULL){
  if (length (freqs) != length (amps)) stop ('Must specify same number of frequencies and amplitudes.')
  if (length (freqs) != length (phases)) stop ('Must specify same number of frequencies and initial phases.')
  amps = abs(amps)
  t = seq (0, dur/1000, 1/fs)  
  n = length (freqs)
  waves = matrix (0, length (t), n)

  if (is.null(colors)) colors = 1:n
  if (length (colors) < n) colors = rep (colors, 100)
  
  for (i in 1:n) waves[,i] = amps[i] * sin (2 * pi * t * freqs[i] + phases[i]) 
  
  if (sum == TRUE) waves = cbind (waves, rowSums (waves))
  
  if (show == TRUE){
    oldpar = par()
    if (sum == TRUE) par (mfrow = c(2,1))
    plot (t*1000,waves[,1], type = 'l', ylab = 'Amplitude', xlab = 'Time (ms)', ylim = c(-max(amps),max(amps)), lwd = 2,xaxs='i')
    abline (h = 0)
    if (n > 1) for (i in 2:n) lines (t*1000, waves[,i], type = 'l', col = colors[i], lwd = 2)
    if (sum == TRUE){
      plot (t*1000,waves[,1+n], col = 1, lwd = 2, type = 'l',xaxs='i', ylab = 'Amplitude', xlab = 'Time (ms)', 
      ylim = range (waves[,1+n]))
      abline (h = 0)
    }
	suppressWarnings (par (oldpar))
  }
  waves = data.frame (time = t*1000, waves)
  colnames(waves)[2:(n+1)] = paste (rep('wave',n), 1:n, sep='')
  if (sum) colnames(waves)[ncol(waves)] = 'sum'
  invisible (waves)
}


sinusoids = function (freqs, amps = rep(1, length(freqs)), dur = 50, phases = rep(0, length(freqs)), fs = 10000, sum = FALSE, show = FALSE, colors = NULL){
  cl = match.call()
  args = sapply (2:length(cl), function(x) cl[[x]])
  names(args) = names(cl)[-1]
  do.call (sinusoid, args)
}

