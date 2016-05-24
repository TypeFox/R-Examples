WiSEConfidenceRegion <-
function(X, Y, plot=TRUE, ...){
  ##Check the Y and X are from WiSE.bootstrap and took only a vector of data
  if(inherits(Y, "WiSEBoot")==FALSE){
    stop("Y should be of class wise_bootstrap (i.e. generated from the WiSE.bootstrap function)")
  }else if(dim(Y$DataWavelet)[2]>1){
    stop("Y should have data input (X) as a vector.")
  }

  if(inherits(X, "WiSEBoot")==FALSE){
    stop("X should be of class wise_bootstrap (i.e. generated from the WiSE.bootstrap function)")
  }else if(dim(X$DataWavelet)[2]>1){
    stop("X should have data input (X) as a vector.")
  }

  ##Y and X Need the same selected J0
  if(Y$MSECriteria[which.min(Y$MSECriteria[,2]), 1] != X$MSECriteria[which.min(X$MSECriteria[,2]), 1] ){
    stop("X and Y should have the same selected wavelet smooth level.")
  }
  ##Need to be of the same data length
  if(dim(Y$BootWavelet)[2]!=dim(X$BootWavelet)[2]) {
    stop("X and Y should have input vectors of the same lengths.")
  }
  
  ##Need to have the same number of bootstrap samples
  if(dim(Y$BootWavelet)[1]!=dim(X$BootWavelet)[1]) {
    stop("X and Y should have the same number of bootstrap samples.")
  }
  
  

  ##Check plot##
  if(!(plot %in% c(TRUE, FALSE))){
    stop("plot is logical.")
  }
  dots <- list(...)


  ##Parameters/quantities from the original bootstrap simulation
  J0PlusOne <- Y$MSECriteria[1,which.min(Y$MSECriteria[,2])] #selected smooth level + 1
  estYWavelet <- Y$DataWavelet[seq(2:2^(J0PlusOne + 1)), 1]       #keep coarse coefficients for Y data
  estXWavelet <- X$DataWavelet[seq(2:2^(J0PlusOne + 1)), 1]       #keep coarse coefficients for X data
  bootEstYWavelet <- Y$BootWavelet[ , seq(2:2^(J0PlusOne + 1)), 1]       #keep coarse coefficients for Y boot samples
  bootEstXWavelet <- X$BootWavelet[ , seq(2:2^(J0PlusOne + 1)), 1]       #keep coarse coefficients for X boot samples
  R <- dim(bootEstYWavelet)[1]


  ##Calculate parameter estimates for the data
  meanYWavelet <- sum(estYWavelet)/length(estYWavelet)  #mean of the coarse coefficients for Y data
  meanXWavelet <- sum(estXWavelet)/length(estXWavelet)  #mean of the coarse coefficients for X data

  dataSlope <- sum((estYWavelet - meanYWavelet)*(estXWavelet - meanXWavelet))/sum((estXWavelet - meanXWavelet)^2)
  dataIntercept <- meanYWavelet - dataSlope*meanXWavelet


  ##Calculate bootstrap parameter estimates
  bootSlope <- rep(NA, R)
  bootIntercept <- rep(NA, R)
  for(r in 1:R){
    meanBootYWavelet <- sum(bootEstYWavelet[r, ])/length(bootEstYWavelet[r, ])
    meanBootXWavelet <- sum(bootEstXWavelet[r, ])/length(bootEstXWavelet[r, ])
    bootSlope[r] <- (sum((bootEstYWavelet[r, ]-meanBootYWavelet)*(bootEstXWavelet[r, ]-meanBootXWavelet))/
                     sum((bootEstXWavelet[r, ]-meanBootXWavelet)^2) )
    bootIntercept[r] <- meanBootYWavelet - bootSlope[r]*meanBootXWavelet
  }


  if(plot==TRUE){
    if(is.null(dots$main)==TRUE){dots$main=paste("Bootstrap Parameter Estimates with Selected Threshold Level", J0PlusOne-1)}
    if(is.null(dots$xlab)==TRUE){dots$xlab=expression(paste("Intercept (", alpha, ")", sep=''))}
    if(is.null(dots$ylab)==TRUE){dots$ylab=expression(paste("Slope (", beta, ")", sep=''))}
    if(is.null(dots$xlim)==TRUE){dots$xlim=c(min(bootIntercept, dataIntercept), max(bootIntercept, dataIntercept))}
    if(is.null(dots$ylim)==TRUE){dots$ylim=c(min(bootSlope, dataSlope), max(bootSlope, dataSlope))}
    if(is.null(dots$pch)==TRUE){dots$pch=20}
    dots$x <- bootIntercept
    dots$y <- bootSlope
    do.call("plot", dots)
    points(dataIntercept, dataSlope, pch=20, cex=2, col="red")
  }


  structure(invisible(list(dataSlope=dataSlope, dataIntercept=dataIntercept, 
                           bootSlope=bootSlope, bootIntercept=bootIntercept)), class="WiSEConfidenceRegion")
}
