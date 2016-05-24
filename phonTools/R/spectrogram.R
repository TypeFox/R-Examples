# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

spectrogram = function (sound, fs = 22050, windowlength = 5, timestep = -1000,
padding = 10, preemphasisf = 50, maxfreq = 5000, colors = TRUE, 
dynamicrange = 50, nlevels = dynamicrange, maintitle = "", show = TRUE, 
window = 'kaiser', windowparameter = 3, quality = FALSE){

  if (class(sound) == "ts") fs = frequency(sound)
  if (class(sound) == "sound"){
    fs = sound$fs
    sound = sound$sound
  }
   
  n = ceiling((fs/1000) * windowlength)     
  if (n%%2) n = n + 1

  if (timestep > 0) timestep = floor(timestep/1000 * fs)
  if (timestep <= 0) timestep = floor (length(sound) / -timestep)
  if (preemphasisf > 0) sound = preemphasis (sound, preemphasisf, fs)

  #sound = c(rep(0, floor(n / 2)), sound, rep(0, floor(n / 2)))
  spots = seq (floor(n / 2), length(sound)-n, timestep)
    
  padding = n*padding
  if ((n + padding)%%2) padding = padding + 1
  N = n + padding

  spect = sapply (spots,function(x){
    tmp = sound[x:(x+n-1)] * windowfunc(sound[x:(x+n-1)], window, windowparameter);
    tmp = c(tmp, rep(0, padding));
    tmp = tmp - mean(tmp);
    tmp = fft (tmp)[1:(N/2+1)]; 
    tmp = abs(tmp)^2;
    tmp = log(tmp, 10) * 10;
  })
  spect = t(spect)
  for (i in 1:nrow(spect)) spect[i,1] = min(spect[i,-1])
	
  hz = (0:(N/2)) * (fs/N)
  times = spots * (1000/fs)
  rownames(spect) = as.numeric (round(times, 2))
  colnames(spect) = as.numeric (round(hz, 2))

  if (colors == 'alternate') colors = c('black','red','orange','yellow','white')
  if (maxfreq > (fs/2)) maxfreq = fs/2
  spect = spect - max(spect)

  specobject = list(spectrogram = spect, fs = fs, windowlength = windowlength, 
               timestep = timestep, dynamicrange = dynamicrange, colors = colors, maxfreq=maxfreq)
  class(specobject) = "spectrogram"

  if (show == TRUE) plot(specobject, ylim = c(0, maxfreq), quality = quality)
  invisible (specobject)
} 


