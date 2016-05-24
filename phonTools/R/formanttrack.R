# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


formanttrack = function (sound, timestep = 5, windowlength = 30, 
formants = 5, cutoff = 5000, minformant = 200, maxbw = 600, 
fs = 22050, show = TRUE, periodicity = .5, returnbw = FALSE){

  if (class(sound) == "ts") fs = frequency(sound)
  if (class(sound) == "sound") {
    fs = sound$fs
    sound = sound$sound
  }     
  if (!is.numeric(sound)) stop ('Sound must be numeric.')
  if (timestep < 0) stop ('Timestep must be positive.')
  if (windowlength < timestep) stop ('Window length must be greater than or equal to timestep.')
  if (cutoff > fs/2) stop ('Cutoff frequency must be less than or equal to the Nyquist frequency.')

  if (cutoff != fs/2){
    sound = resample (sound, cutoff*2, fs, synthfilter = TRUE)
    fs = cutoff*2
  }
  T = 1 / fs
  stepsize = round ((timestep / 1000) / T)
  half = ceiling (windowlength/1000 * fs)/ 2
  spots = seq (1+half, length(sound)-half, stepsize)

  ffs = NULL
  bws = NULL
  total = 1
  for (i in 1:(length(spots)-1)){
    section = sound[(spots[i]-half):(spots[i]+half)]
    tmpacf = pitchtrack (section, timestep=0, fs=fs)[2] 
    if (tmpacf>=periodicity){
      tmp = findformants (section, fs = cutoff*2, maxbw = maxbw, minformant = minformant, 
      verify = F, coeffs = formants*2+3)[1:formants,]
      tmp[is.na(tmp[,1]),] = 0
      ffs = rbind (ffs, c(0,tmp[,1]))
      bws = rbind (bws, tmp[,2])
      ffs[total,1] = spots[i] * (1000/fs)
      total = total + 1
    }
  }
  colnames(ffs) = c('time',paste ('f', 1:formants, sep=''))
  rownames(ffs) = 1:nrow(ffs)
  
  if (returnbw){
    ffs = cbind (ffs, bws)
    colnames(ffs) = c('time',paste ('f', 1:formants, sep=''),paste ('b', 1:formants, sep=''))
  }
    tmp = data.frame (ffs)
  
  if (show == TRUE){
    spectrogram (sound, colors = FALSE, quality = FALSE, fs = fs)
    for (i in 2:(formants+1)) points (ffs[ffs[,i]!=0,1], ffs[ffs[,i]!=0,i], pch = 16, col = i-1, cex = .75)
  }
  invisible (tmp)
}

