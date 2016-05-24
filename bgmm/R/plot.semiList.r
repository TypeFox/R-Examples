plotGIC <- function(models, penalty=2, plot.it=TRUE, ...) {
  if (!("mModelList" %in% class(models))) {
    if (plot.it) {
      dotchart(models, xlab=paste("GIC"),pch=19)
#      abline(v = min(models), col='grey', lty=2)
    }
    return(models)
  }
  if (length(models$models)==length(models$kList) || length(models$kList)==1) {
     # only one variable
    values = sapply(models$models, getGIC, p=penalty, ...)
    tnames = models[[3]]
    names(values) = tnames
    if (plot.it) {
      ind = which.min(values)
      tnames[ind] = paste("*",tnames[ind])
      tnames[-ind] = paste("  ",tnames[-ind])
      dotchart(values, labels=tnames, xlab=paste("GIC, penalty =",penalty),pch=19)
#      abline(v = min(values), col='grey', lty=2)
    }
    return(values)
  } else {
    values = sapply(models$models, getGIC, p=penalty, ...)
    tmp= matrix(values, ncol=length(models$kList))
    rownames(tmp) = sapply(strsplit((models[[3]])[1:(length(models[[3]])/length(models$kList))]," "), tail,1)
    colnames(tmp) = paste("k=",models$kList,sep="")
    if (plot.it) {
      dotchart(tmp,  xlab=paste("GIC, penalty =",penalty),pch=19)
      mm = (nrow(tmp)+2)*(ncol(tmp) - which.min(apply(tmp,2,min))) + which.min(apply(tmp,1,min))
      mtext(side=2,at=mm,"*",line=1)
#      abline(v = min(tmp), col='grey', lty=2)
    }
    return(tmp)
  }
}



plot.mModelList <- function(x, ...) {
  d = x$models[[1]]$d
  X = x$models[[1]]$X
  knowns = x$models[[1]]$knowns
  B = x$models[[1]]$B
#  if (d > 2) 
#    stop("PLOT SUPPORTS ONLY 1D and 2D data")
  if (d==1) 
    plotList.1d(X, knowns, B, x, ...)
  if (d > 1) 
    plotList.2d (X, knowns, B, x, ...)
}

plotList.2d <-function(X, knowns, B, models2d, ...) {
  if (length(models2d$models)==length(models2d$kList) || length(models2d$kList)==1) {
     # only one variable
    gridwi = floor(sqrt(length(models2d[[2]])))
    gridhe = ceiling(length(models2d[[2]])/gridwi)
  } else {
    gridwi = length(models2d$kList)
    gridhe = length(models2d$models)/length(models2d$kList)
  }
  par(mfrow=c(gridwi, gridhe), mar=c(2,2,2,1))
  for (i in seq_along(models2d[[2]])) 
    plot.2d(X, knowns, map(B), models2d[[1]][[i]], main=models2d[[3]][[i]], ...) 
}

plotList.1d <-function(X, knowns, B, models1d, ...) {
  if (length(models1d$models)==length(models1d$kList) || length(models1d$kList)==1) {
     # only one variable
    gridwi = floor(sqrt(length(models1d[[2]])))
    gridhe = ceiling(length(models1d[[2]])/gridwi)
  } else {
    gridwi = length(models1d$kList)
    gridhe = length(models1d$models)/length(models1d$kList)
  }

  par(mfrow=c(gridwi, gridhe), mar=c(2,2,2,1))
  for (i in seq_along(models1d[[2]])) 
    plot.1d(X, knowns, map(B), models1d[[1]][[i]], main=models1d[[3]][[i]], ...) 
}

