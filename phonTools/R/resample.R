# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


resample = function (sound, newfs, oldfs, precision = 50, filterorder = 200, synthfilter = FALSE){
  soundout = 0; tsout = 0;
  if (class(sound) == "ts"){
    fs = frequency(sound)
    tsout = 1
  } 
  if (class(sound) == "sound") {
    soundout = 1
    oldsound = sound
    oldfs = sound$fs
    sound = sound$sound
  } 
  ratio = oldfs / newfs 
  if (ratio > 1 & !synthfilter) sound = FIRfilter (sound, to = newfs/2, fs = oldfs, order = filterorder)
  if (ratio > 1 & synthfilter) sound = synthfilter (sound, band = c(0,newfs/2), fs = oldfs)
  
  newtime = seq (1, length(sound)+1, by = ratio)   
  nearest = round (newtime)                              
  offset = newtime - nearest                                 
  
  sound = c(rep(0,precision), sound, rep(0,precision+1))
  y = newtime * 0
  
  for (i in -precision:precision)
    y = y + sound[nearest+precision+i] * sinc(offset - i, normalized = TRUE)
  
  if (ratio < 1 & !synthfilter) y = FIRfilter (y, to = oldfs/2, fs = newfs, order = filterorder)
  if (ratio < 1 & synthfilter) y = synthfilter (y, band = c(0,oldfs/2), fs = newfs)
  
  sound = y / (max(y) * 1.05) 
  if (soundout == 1)  sound = makesound (sound, filename = oldsound$filename, fs = newfs)
  if (tsout == 1)  sound = ts (sound, frequency = newfs, start = 0)
  return (sound)   
}

