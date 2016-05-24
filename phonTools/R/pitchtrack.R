# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


pitchtrack = function (sound, f0range = c(60,400), timestep = 2, fs = 22050, minacf = .5,
                       correction = TRUE, show = TRUE, windowlength = 50, addtospect = FALSE){
  if (length(f0range) != 2) stop ('A highest and lowest f0 must be specified.')
  if (f0range[2] < f0range[1]) stop ('Maximum f0 must be greater than minimum f0.')
  if (minacf < 0) stop ('minacf must be positive.')
  if (timestep<0) stop ('Timestep must be positive.')
 
  if (class(sound) == "ts") fs = frequency(sound)
  if (class(sound) == "sound") {
    fs = sound$fs
    sound = sound$sound
  }  
  if (!is.numeric(sound)) stop ('Sound must be numeric.')

  T = 1 / fs
  stepsize = round ((timestep / 1000) / T)
  minlag = round (1 / (T * f0range[2]))
  maxlag = round (1 / (T * f0range[1] ))

  if (timestep>0){
    half = round( ceiling (windowlength/1000 * fs)/ 2)
    spots = seq (1+half, length(sound)-half, stepsize)
    corr = rep (0, length(spots))
    lag = rep (0, length(spots))
    for (i in 1:length(spots)){
      section = sound[(spots[i]-half):(spots[i]+half)] 
      acf = fastacf (section, lag.max = maxlag, show = F, correct = correction)    
      peaks = peakfind (acf[,2], show = FALSE)
      lag[i] = peaks[order(acf$acf[peaks], decreasing = TRUE)[1]]
      if (is.na(lag[i])) lag[i] = 0
      if (lag[i]!=0) corr[i] = acf$acf[lag[i]]
      if (lag[i]==0) corr[i] = 0
    }
    lag[lag < minlag | lag > maxlag] = 0
    spots = spots * T * 1000
    f0 = 1 / (lag * T)
    f0 [f0 == Inf] = 0
    f0 [f0 == fs] = 0

    use = (f0 != 0 & corr > minacf)
    spots = spots[use]; corr = corr[use]; f0 = f0[use];

    if (show & !addtospect) plot (spots, f0, cex = 1.2*corr, pch = 16, col = 4, ylab = 'f0 (Hz)', 
    xlab = 'Time (ms)', ylim = c(0,f0range[2]))
  
    if (addtospect) points (spots, f0*10, col = 6, pch = 16, cex = 1.2*corr)
    output =  data.frame (time = round(spots,1), f0 = round(f0,2), acf = round(corr,4))
  }
  if (timestep==0){
    acf = fastacf (sound, lag.max = maxlag, show = F, correct = correction)    
    peaks = peakfind (acf[,2], show = FALSE)
    lag = peaks[order(acf$acf[peaks], decreasing = TRUE)[1]]
    if (is.na(lag)) lag = 0
    if (lag < minlag) lag = 0
    f0 = 1 / (lag * T)
    f0 [f0 == Inf] = 0
    acf = acf$acf[lag]
    if (lag==0) acf = 0
    output = (c(f0, acf))
	names (output) = c('f0','acf')
  } 
  invisible (output) 
}
  
  