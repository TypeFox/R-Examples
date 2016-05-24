# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


imputeformants = function (ffs, speaker, vowel){
  if (!is.numeric (ffs)) stop ('Non-numeric formant frquency values provided.') 
  
  vowel = as.factor (vowel)
  speaker = as.factor (speaker)
  lvs = levels (vowel)
  lsp = levels (speaker)

  nvs = length (lvs)
  nsp = length (lsp)

  if (max(ffs) > 20) ffs = log (ffs)
  omit = !(ffs == -Inf)
  mod = lm (ffs[omit] ~ vowel[omit] + speaker[omit])
  cffs = mod$coefficients
  
  vcffs = as.numeric (c(cffs[1], cffs[1]+cffs[2:nvs]))
  pcffs = c(0,cffs[(nvs+1):(nvs+nsp-1)])

  ffs = exp(ffs) 
  for (i in 1:length (ffs))
    if (ffs[i] == 0) 
      ffs[i] = exp (vcffs[lvs == vowel[i]] + pcffs[lsp == speaker[i]])
  
  return (round(ffs))
}

