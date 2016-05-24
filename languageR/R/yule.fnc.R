`yule.fnc` <-
function(spect) {
  N = sum(spect$frequency * spect$freqOfFreq)
  return(10000 * (sum(spect$frequency^2 * spect$freqOfFreq) - N)/N^2)
}

