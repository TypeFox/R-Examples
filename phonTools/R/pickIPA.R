# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


pickIPA = function (vowels, n = 0, xsampa = FALSE, description = FALSE, verify = TRUE){
  if (missing(vowels) & n == 0) stop ('Please provide either a vowels vector or the number of vowels to select.')
  vector = FALSE
  if (!missing(vowels)) vector = TRUE
  if (vector){
    vowels = as.factor (vowels)
    vtypes = levels(vowels)
    nvowels = ntypes (vowels)
    n = nvowels
  }
  IPA = ipainfo()[c(1,4,2,5)]
  plot (IPA[[2]]$frontness+(IPA[[2]]$rounded*.25), IPA[[2]]$height, pch = IPA[[1]],cex = 2, 
  col = 1+IPA[[2]]$rounded*3, xlab = '',ylab = '', xaxt = 'n', yaxt = 'n',xlim = c(.75,3.5),
  ylim=c(0.75,4.25))
  axis (side = 1, at = c(1.15,2.15,3.15), c('Front','Mid','Back'), cex.axis = 1.3)
  axis (side = 2, at = c(1,2,3,4), c('open','open-mid','close-mid','close'), cex.axis = 1.3)

  selected = rep(0,n)
  for (i in 1:n){ 
    if (vector) cat ('Please select ->  ', vtypes[i], '\n\n')
    if (!vector) cat ('Please select vowel ', i, '\n\n')
    flush.console()
    selected[i] = identify (IPA[[2]]$frontness+(IPA[[2]]$rounded*.25), IPA[[2]]$height,'', n = 1)  
  }

  if (verify == TRUE) plot (1:n, rep(1,n), pch = IPA[[1]][selected], ylab='',yaxt='n', xlab='Selection',cex = 2)
  if (vector == TRUE) selected = selected[as.numeric(vowels)]
  out = list (IPA = IPA[[1]][selected])
  if (xsampa == TRUE) out$xsampa = IPA[[4]][selected]
  if (description == TRUE) out$description = IPA[[3]][selected]
  out
}
