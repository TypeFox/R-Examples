# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


powertrack = function (sound, timestep = 5, windowlength = 30, 
                       fs = 22050, show = TRUE, zeromax = TRUE, ...){
  if (class(sound) == "ts") fs = frequency(sound)
  if (class(sound) == "sound") {
    fs = sound$fs
    sound = sound$sound
  }  
  if (!is.numeric(sound)) stop ('Sound must be numeric.')
  if (timestep < 0) stop ('Timestep must be positive.')
  if (windowlength < timestep) stop ('Window length must be greater than or equal to timestep.')

  T = 1 / fs
  stepsize = round ((timestep / 1000) / T)
  
  half = round(ceiling (windowlength/1000 * 22050)/ 2)
  
  spots = seq (1+half, length(sound)-half, stepsize)
  power = rep (0, length(spots))
  
  for (i in 1:length(spots)){
    section = sound[(spots[i]-half):(spots[i]+half)] 
    power[i] = mean ((section*windowfunc(section))^2)    
  }
  use = (power != 0)
  spots = spots[use]
  power = power[use]
  
  power = log(power, 10)*10
  if (zeromax) power = power - max(power)
  
  time = spots * (1000/fs)
  tmp = data.frame (time = time, power = power)
  
  if (show == TRUE) plot(tmp$time, tmp$power, xlab = 'Time (ms)', ylab = 'Power (dB)', 
                    type = 'l', ylim = c(min(power)-1, 2), lwd = 2, col = 4, ...) 
  invisible (tmp)
}

