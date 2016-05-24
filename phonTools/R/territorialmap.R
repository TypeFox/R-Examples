# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

territorialmap = function (template, show = TRUE){
  if (class(template)!='template') stop ('Templates must be created with createtemplate().')
  means = template$means
  covariance = template$covariance
  xlim = template$ranges[1,]
  ylim = template$ranges[2,]
  tempobject = TRUE
  pols = list()

  for (j in 1:nrow(means)){
    coeffs = matrix (0, nrow(means),2)
    plot (means[,1:2], pch = rownames (means), cex = 2, 
          ylim = ylim, xlim = xlim, main = rownames (means)[j])
    for (i in (1:nrow(means))[-j]){ 
      coeffs[i,] = ldboundary (means[j,1:2],means[i,1:2],covariance[1:2,1:2], add = T)
    }
    coeffs = coeffs[-j,]
    coeffs = rbind (coeffs, c(ylim[1],0), c(ylim[2],0))
    b = (ylim[2] - ylim[1]) / .001
    a1 = ylim[1] + b*xlim[1]
    a2 = ylim[2] + b*xlim[2]
    coeffs = rbind (coeffs, c(a1, -b))
    coeffs = rbind (coeffs, c(a2, -b))
    for (i in nrow(coeffs):(nrow(coeffs)-3)) abline (coeffs[i,])

    x = NULL; y = NULL;
    count = 0
    for (i in 1:(nrow(coeffs)-1)){
      for (k in (i+1):nrow(coeffs)){
        count = count + 1
        x = c(x, (coeffs[i,1] - coeffs[k,1]) / (coeffs[k,2] - coeffs[i,2]))
        y = c(y, coeffs[i,1] + coeffs[i,2] * x[count])
        points (x[count],y[count],cex=2, col=4, pch=15 + count%%4)
      }
    }
    selections = NULL
    tmp = NULL
    count = 0
    while (sum (tmp == selections) < 2 & !is.null(dev.list())){
      count = count + 1
      tmp = identify (x,y,count,n = 1, col = 2)
      points (x[tmp], y[tmp], col = 2, pch=16, cex=2)
      selections = c(selections, tmp)
    }
    pols[[j]] = cbind (x[selections],y[selections])
  }
  if (is.null(dev.list())) stop ('Plot window closed.')
  template$territory = pols
  if (show) plot (template, territorial = TRUE)
  return (template)
} 


