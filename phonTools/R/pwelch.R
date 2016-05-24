# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

pwelch = function (sound, points = 0, overlap = 0, padding = 0, window = 'hamming', 
                   show = TRUE, fs = 1, preemphasisf = 0, zeromax = TRUE, type,...){
  if (class(sound) == "sound") {
    fs = sound$fs
    sound = sound$sound
  }
  if (class(sound) == "ts") fs = frequency(sound)

  if (preemphasisf > 0) sound = preemphasis (sound, preemphasisf, fs)
  
  n = length (sound)
  if (points == 0) points = ceiling (n / 10)
  
  sound = c(sound, rep(0, points))
  
  spots = seq (1, n, points - overlap)
  
  if ((points+padding) %% 2 == 1) padding = padding + 1
  n = points + padding
  
  magnitude = rep (0, n)
  for (i in 1:length(spots)){
    tmp = sound[spots[i]:(spots[i]+points-1)] * windowfunc(points, type = window)
    tmp = c(tmp, rep (0, padding))
    tmp = fft(tmp)
    tmp = tmp * Conj (tmp)
    magnitude = magnitude + tmp
  }
  magnitude = magnitude / length(spots)
  magnitude = magnitude[1:(n/2+1)]
  magnitude = abs(magnitude)
  dB = log (magnitude, 10) * 20
  if (zeromax == TRUE) dB = dB - max (dB)
  
  if (fs > 1) hz = seq (0, fs/2, length.out = (n/2)+1)
  if (fs == 1) hz = seq (0, .5, length.out = (n/2)+1)
  
  if (missing(type)) type = "l"
  
  if (fs > 1) xlab = 'Frequency (Hz.)'
  if (fs == 1) xlab = 'Frequency / Sampling Freq.'

  if (show == TRUE) plot (hz, dB, type = type, ylab = 'Power (dB.)', xlab = xlab, xaxs = 'i' ,...)  
  
  invisible (cbind (hz, dB))
}

