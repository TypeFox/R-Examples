`zipf.fnc` <-
function(text, plot = FALSE) {
  spectrum = spectrum.fnc(text)
  spectrum = spectrum[nrow(spectrum):1,]
  spectrum$rank = cumsum(spectrum$freqOfFreq)
  if (plot) 
    plot(log(spectrum$rank), log(spectrum$frequency), 
      xlab = "log rank", ylab = "log frequency", 
      type = "S")
  return(spectrum)
}

