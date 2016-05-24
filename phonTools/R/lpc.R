# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


lpc = function (sound, order = round(fs/1000) + 3, fs = 10000, show = FALSE, add = FALSE, preemph = TRUE){
  if (class(sound) == "ts") fs = frequency(sound)
  if (class(sound) == "sound"){
    fs = sound$fs
    sound = sound$sound
  }

  if (!is.numeric(sound)) stop("Input must be numeric.")
  if (preemph == TRUE) sound = preemphasis(sound, fs = fs)
  n = length (sound)              
  sound = sound - mean(sound)
  
  sound = sound * windowfunc (sound)
  sound = c(sound, rep(0, order))              
  
  predictors = t(sapply(seq(1, n, 1), function(x) sound[(x):(x + order)]))  
  y = sound[1:n]           
  r = y %*% predictors              
  
  tmp = c(rev(r), r[-1])
  w = t(sapply(seq(order+1, 2, -1), function(x) tmp[(x):(x+order-1)]))       

  coeffs = -r[2:(order+1)] %*% solve(w)
  coeffs = c(1, coeffs)
  
  if (show == TRUE & add == TRUE) 
    freqresponse(1, coeffs, fs = fs, add = add)
  if (show == TRUE & add == FALSE) {
    spectralslice(sound, fs = fs, col = 4)
    freqresponse(1, coeffs, fs = fs, add = TRUE)
  }
  coeffs
}              

