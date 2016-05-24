# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


xsampatoIPA = function (vowels, chart = FALSE, verify = FALSE){
  oldpar = par()
  if (missing(vowels) | chart == TRUE){ 
    par (mfrow = c(1,2))
    IPA = ipainfo()[c(1,4,5)]
    plot (IPA[[2]]$frontness+(IPA[[2]]$rounded*.25), IPA[[2]]$height, pch = IPA[[1]],cex = 1.5, 
          col = 1+IPA[[2]]$rounded*3, xlab = '',ylab = '', xaxt = 'n', yaxt = 'n',xlim = c(.75,3.5),ylim=c(0.75,4.25))
    axis (side = 1, at = c(1.15,2.15,3.15), c('Front','Mid','Back'))
    axis (side = 2, at = c(1,2,3,4), c('open','open-mid','close-mid','close'))
  
    plot (IPA[[2]]$frontness+(IPA[[2]]$rounded*.25), IPA[[2]]$height, type = 'n', cex = 1.5, 
          col = 1+IPA[[2]]$rounded*3, xlab = '',ylab = '', xaxt = 'n', yaxt = 'n',xlim = c(.75,3.5),ylim=c(0.75,4.25))
    text (IPA[[2]]$frontness+(IPA[[2]]$rounded*.25), IPA[[2]]$height, label = IPA[[3]],cex = 1.25, 
          col = 1+IPA[[2]]$rounded*3)
    axis (side = 1, at = c(1.15,2.15,3.15), c('Front','Mid','Back'))
    axis (side = 2, at = c(1,2,3,4), c('open','open-mid','close-mid','close'))
	suppressWarnings (par (oldpar))
  }
  if (!missing(vowels)){ 
    IPA = ipainfo()[c(1,5)]
    out = as.hexmode(rep (0, length (vowels)))
    for (i in 1:length(IPA[[2]]))
      if (sum(IPA[[2]][i] == vowels)>0) out[IPA[[2]][i] == vowels] = IPA[[1]][i]
    if (verify == TRUE){
	  par (mfrow = c(1,1))
      plot (1:length(out), rep(1,length(out)), pch = out, ylab='',yaxt='n', xlab='Selection',cex = 2,xaxt = 'n')
	  axis (side = 1, at  = 1:length(out))
      suppressWarnings (par (oldpar))
    }
    return(out)
  }
 
}

