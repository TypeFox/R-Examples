# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


playsound = function (sound, path = 'default', fs = 10000, erase = TRUE){
  if (class (sound) != 'sound' & class (sound) != 'ts'){
    sound = sound / (max(sound)*1.05)
    sound = makesound (sound, 'play_tmp.wav', fs = fs)
  }
  writesound (sound, 'play_tmp.wav')

  pathout = FALSE
  if (path == 'default') path = "C:\\Program Files (x86)\\VideoLAN\\VLC\\vlc.exe"
  if (path == 'pick'){
    path = file.choose()
    pathout = TRUE
  }
  
  if (!file.exists (path)) stop ('VLC player is not installed or incorrect path.')
  system(paste(shQuote(path), '--play-and-exit','play_tmp.wav'),ignore.stderr=TRUE)
  if (erase) unlink ('play_tmp.wav')
  if (pathout) return (path)
}

  
play = function (sound, path = 'default', fs = 10000, erase = TRUE){
  cl = match.call()
  args = sapply (2:length(cl), function(x) cl[[x]])
  names(args) = names(cl)[-1]
  do.call (playsound, args)
}

