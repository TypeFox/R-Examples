# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


normalize = function (formants, speakers, vowels, method = 'neareyE', corners = NULL){
  if (is.null(ncol(formants))) stop("At least two formants must be provided (i.e. F1, F2, ...)")
  if (length(speakers) != nrow (formants)) stop('Speaker vector length does not match formant data length.')
  if (length(vowels) != nrow (formants)) stop('Formant vector length does not match formant data length.')
  if (!(method %in% c('barreda','neareyE', 'neareyI','lobanov','wandf'))) stop ('Invalid method selected. See help file for available methods.')
  
  speakers = as.factor (as.character(speakers))
  vowels = as.factor (as.character(vowels))
  speakersf = levels (speakers)
  vowelsf = levels (vowels)
  if (method == 'wandf') 
  if (sum (corners %in% vowelsf) != 2 | is.null(corners)) stop ('Please provide two corner vowels which are present in the vowel vector.')

  ns = length (speakersf) 
  nv = length (vowelsf)
  nffs = ncol (formants)  
  
  if (method == 'neareyE'){
    if (max(formants) > 30) formants = log (formants)
    meanffs = rowSums(formants)/ncol(formants)
    for (j in 1:ns){
      temp = (speakers == speakersf[j])
      psi = mean (tapply(meanffs[temp], vowels[temp], mean))
      formants[temp, ] = formants[temp, ] - psi
  }}
  if (method == 'neareyI'){
    if (max(formants) > 30) formants = log (formants)
    for (j in 1:ns){
      temp = (speakers == speakersf[j])
      mff = NULL
      for (i in 1:nffs) mff = c(mff, mean (tapply(formants[temp,i], vowels[temp], mean)))
      mffs = matrix (mff, nrow (formants[temp,]), nffs, byrow = TRUE)
      formants[temp, ] = formants[temp, ] - mffs
  }}
  if (method == 'lobanov'){
    for (j in 1:ns){
      temp = (speakers == speakersf[j])
      mff = NULL
      sdff = NULL
      for (i in 1:nffs){ 
         mff = c(mff, mean (tapply(formants[temp,i], vowels[temp], mean)))
         sdff = c(sdff, sd(tapply(formants[temp,i], vowels[temp], mean)))
      }
      mffs = matrix (mff, nrow (formants[temp,]), nffs, byrow = TRUE)
      sdffs = matrix (sdff, nrow (formants[temp,]), nffs, byrow = TRUE)
      formants[temp, ] = (formants[temp, ] - mffs) / sdffs
  }}
  if (method == 'barreda'){
    if (max(formants) > 30) formants = log (formants)
    for (j in 1:ns){
      temp = (speakers == speakersf[j])
      mff = NULL
      sdff = NULL
      for (i in 1:nffs){ 
         mff = c(mff, sum(range (tapply(formants[temp,i], vowels[temp], mean)))/2)
         sdff = c(sdff, mean((mff[i] - tapply(formants[temp,i], vowels[temp], mean))^2)^.5)
      }
      mffs = matrix (mff, nrow (formants[temp,]), nffs, byrow = TRUE)
      sdffs = matrix (sdff, nrow (formants[temp,]), nffs, byrow = TRUE)
      formants[temp, ] = (formants[temp, ] - mffs) / sdffs
  }}  
  if (method == 'wandf'){
    formants = formants[,1:2]
    for (j in 1:ns){
      temp = (speakers == speakersf[j])
      iyf1 = mean (formants[temp & vowels == corners[1], 1])
      iyf2 = mean (formants[temp & vowels == corners[1], 2])
      ahf1 = mean (formants[temp & vowels == corners[2], 1])
      ahf2 = mean (formants[temp & vowels == corners[2], 2])
     
      sf1 = (iyf1 + ahf1 + iyf1) / 3
      sf2 = (iyf2 + ahf2 + iyf1) / 3
     
      formants[temp, 1] = formants[temp, 1] / sf1
      formants[temp, 2] = formants[temp, 2] / sf2
  }}  
  
  output = data.frame (formants, speaker = speakers, vowel = vowels)
  return (output)
}

