# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


findformants = function (sound, fs = 10000, coeffs = NULL, maxbw = 600, 
minformant = 200, verify = TRUE, showbws = FALSE, showrejected = TRUE){
  if (missing (sound)) sound = 1

  if (class(sound) == "ts") fs = frequency(sound)
  if (class(sound) == "sound") {
    fs = sound$fs
    sound = sound$sound
  } 

  if (is.null(coeffs)) coeffs = lpc (sound, fs = fs)
  if (length(coeffs) == 1) coeffs = lpc (sound, fs = fs, order = coeffs)
  
  roots = polyroot (rev(coeffs))
  angs = atan2 (Im(roots), Re(roots))
  formants = round (angs * (fs/(2*pi)), 2)
  nums = order (formants)
  formants = formants[nums]
  bws = -(fs/pi) * log (abs(roots[nums]))
  touse = (bws < maxbw & formants > minformant & formants < fs/2)
  out = data.frame (formant = formants[touse], bandwidth = bws[touse])
  
  if (verify == TRUE){
    multiplot (sizes = c(.7,.3), type = 'c', show = FALSE)
    cols = rep (2:6, 10)
    freqresponse (1, coeffs, fs = fs)
    if (length(sound) > 1) spectralslice (preemphasis(sound,fs=fs), fs = fs, add = TRUE, padding = 0, col = 1, lty = 'dotted')
    for (i in 1:nrow(out)){
      abline (v = out[i,1], lwd = 2, col = cols[i])
      if (showbws == TRUE) abline (v = out[i,1] + out[i,2], lty = 'dotted', col = cols[i])
      if (showbws == TRUE) abline (v = out[i,1] - out[i,2], lty = 'dotted', col = cols[i])
    }    
    if (showrejected == TRUE) abline (v = formants[!touse], lty = 'dotted', lwd = 2)
    plot (roots[nums], xlim = range (-1.1,1.1), ylim = range (-1.1,1.1), pch = 4, lwd = 2, 
          xlab = 'Real', ylab = 'Imaginary', col = !touse)
    sdellipse (means = c(0,0), points = matrix (c(1,0,0,1),2,2), stdev = 1, density = .01)
    abline (h = 0, v = 0, lty = 'dotted')
    tmp = 0
    for (i in 1:length(touse))  
      if (touse[i]){
        tmp = tmp + 1
        points (roots[nums][i], pch = 4, lwd = 2, col = cols[tmp])
      }    
  }
  invisible (out)
}

