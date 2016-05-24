`spectrum.fnc` <-
function(text) {
  tab = table(table(text))
  spectrum = data.frame(frequency = as.numeric(rownames(tab)),
                        freqOfFreq = as.numeric(tab))
  return(spectrum)
}

