# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

mscohere = function (signal1, signal2, points = 0, overlap = 0, padding = 0, 
            pval = TRUE, window = 'hamming', show = TRUE, fs = 1, removeDC = TRUE, 
			ylim, xlim, ind = FALSE, add = FALSE,...){
			
  n1 = length (signal1)
  n2 = length (signal2)
  
  if (n1 != n2) stop ('Signals must be of same length.')
  
  if (points == 0) points = ceiling (n1 / 10)
  
  signal1 = c(signal1, rep(0, points))
  signal2 = c(signal2, rep(0, points))
  
  spots = seq (1, n1, points - overlap)
  
  if ((points + padding) %% 2 == 1) padding = padding + 1
  n = points + padding
  
  XX = rep (0, n)
  YY = rep (0, n)
  XY = rep (0, n)
  
  XXS = NULL
  YYS = NULL
  XYS = NULL
  
  if (removeDC){
    tmp1 = signal1 - mean(signal1)
    tmp2 = signal1 - mean(signal1)
  }
  
  for (i in 1:length(spots)){
    tmp1 = signal1[spots[i]:(spots[i]+points-1)]
    tmp2 = signal2[spots[i]:(spots[i]+points-1)] 
    
    tmp1 = tmp1 * windowfunc(points, type = window)
    tmp2 = tmp2 * windowfunc(points, type = window)
    tmp1 = c(tmp1, rep (0, padding))
    tmp2 = c(tmp2, rep (0, padding))
    
    tmp1 = fft(tmp1) / (1+length(tmp1)-sum(tmp1 == 0))
    tmp2 = fft(tmp2) / (1+length(tmp2)-sum(tmp2 == 0))
    
    if (ind){
      XXS = cbind(XXS, tmp1)
      YYS = cbind(YYS, tmp2)
      XYS = cbind(XYS, tmp1 * Conj(tmp2))
    }
    
    XX = XX + tmp1 * Conj(tmp1)
    YY = YY + tmp2 * Conj(tmp2)
    XY = XY + tmp1 * Conj(tmp2)
  }
  
  XX = XX / length(spots)
  YY = YY / length(spots)
  XY = XY / length(spots)
  
  XX = XX[1:(n/2+1)]
  YY = YY[1:(n/2+1)]
  XY = XY[1:(n/2+1)]
  
  coh = abs(XY)**2 / (abs(XX) * abs(YY))
  hz = seq (0, fs/2, length.out = (n/2)+1)
  
  if (missing(xlim)) xlim = c(0,fs/2)
  if (missing(ylim)) ylim = c(0,max(coh[hz < max(xlim) & hz > min(xlim)]))
  
  p = 1 - 0.05^(1/(length (spots)-1));
  
  tmp = mean (coh[hz < xlim[2] & hz > xlim[1]] > p)
  #if (pval) cat ('% significant coefficients:', tmp, '\n\n')
  
  if (show == TRUE){
    if (fs > 1 & !add) plot (hz, coh, type = 'l', ylab = 'Coherence', xlab = 'Frequency (Hz.)',xlim = xlim, ylim = ylim,xaxs = 'i', yaxs='i',...)  
    if (fs == 1 & !add) plot (hz, coh, type = 'l', ylab = 'Coherence', xlab = 'Frequency / Sampling Freq.',xlim=xlim, ylim = ylim,xaxs = 'i',yaxs='i',...)  
    
    if (fs > 1 & add) lines (hz, coh, ...)  
    if (fs == 1 & add) lines (hz, coh, ...)  
    
    if (pval) abline (h = p, lty = 'dotted', col = 2)
  }
  if (ind) return (list (cbind (hz, coh), list (XXS, YYS, XYS), pval = p))
  if (!ind) invisible (list (data.frame (hz, coh), pval = p))
}


