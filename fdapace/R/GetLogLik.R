# k: input denoting number of components used
# returns -2 times log-likelihood
GetLogLik = function(fpcaObj, k, y = NULL, t = NULL){
  if(fpcaObj$optns$lean == TRUE && (is.null(y) || is.null(t))){
    stop("Option lean is TRUE, need input data y and measurement time list t to calculate log-likelihood.")
  }
  if(fpcaObj$optns$lean == FALSE){ # when input data is in fpcaObj
    y <- fpcaObj$inputData$y
    t <- fpcaObj$inputData$t
  }
  lambda = fpcaObj$lambda[1:k]
  sigma2 = fpcaObj$sigma2
  if(is.null(sigma2) && fpcaObj$optns$dataType == "Dense"){
    ymat = matrix(unlist(y),nrow=length(y), byrow=TRUE)
    sddiag = sqrt(diag(var(ymat)))
    sigma2 = sddiag*1e-4
    sigma2 = ConvertSupport(fromGrid = fpcaObj$obsGrid, toGrid = fpcaObj$workGrid, mu = sigma2)
  }
  logLik = 0
  phi = fpcaObj$phi[,1:k, drop=FALSE]

  if(fpcaObj$optns$dataType %in% c('Dense'
    #, 'DenseWithMV' # need extra imputation step
    )){
  	if(k == 1){
  	  Sigma_y = phi %*% (lambda*diag(k)) %*% t(phi) + sigma2*diag(rep(1,nrow(phi)))
  	} else {
      Sigma_y = phi %*% diag(lambda) %*% t(phi) + sigma2*diag(rep(1,nrow(phi)))
    }
    #detSigma_y = prod(c(lambda,rep(0,nrow(phi)-k))+sigma2)
    detSigma_y = det(Sigma_y)
    if(detSigma_y == 0){
      logLik = NULL
      return(logLik)
    }
    for(i in 1:length(y)){ # TODO: imputation for DenseWithMV needed
      invtempSub = solve(Sigma_y, y[[i]] - fpcaObj$mu)
      logLik = logLik + log(detSigma_y) + invtempSub %*% (y[[i]] - fpcaObj$mu)
    }
    return(logLik)
  } else { # Sparse case
    if(is.null(sigma2)){ sigma2 <- fpcaObj$rho }
    if(fpcaObj$optns$error == TRUE && sigma2 <= fpcaObj$rho){
      # especially for the case when sigma2 is estimated to be <=0 and set to 1e-6
      sigma2 <- fpcaObj$rho
    }
    for(i in 1:length(y)){
      if(length(t[[i]]) == 1){
        phi_i = t(as.matrix(ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = t[[i]],
                               phi = phi)))        
      } else {
        phi_i = ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = t[[i]],
                               phi = phi)
      }
      mu_i = ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = t[[i]],
        mu = fpcaObj$mu)
      if(k == 1){
        Sigma_yi = phi_i %*% (lambda*diag(k)) %*% t(phi_i) + sigma2 * diag(rep(1,length(mu_i)))
      } else{
        Sigma_yi = phi_i %*% diag(lambda) %*% t(phi_i) + sigma2 * diag(rep(1,length(mu_i)))
      }
      detSigma_yi = det(Sigma_yi)
      if(detSigma_yi == 0){
        logLik = NULL
        return(logLik)
      }
      invtempi = solve(Sigma_yi, y[[i]] - mu_i)
      logLik = logLik + log(detSigma_yi) + invtempi %*% (y[[i]] - mu_i)
    }
    return(logLik)
  }
}
