maternCholSolve = function(param, obsCov, coordinates){
  
  cholCovMat = geostatsp::matern(x=coordinates, 
      param=param, type='cholesky')
  detCholCovMat =  attributes(cholCovMat)$logDetHalf

  if(FALSE){ # do this in R
#  covMat = geostatsp::matern(x=coordinates, param=param)	
#  cholCovMat = Matrix::chol(covMat)
    
    
  # cholCovMat %*% t(cholCovMat) = covMat
  
  # cholCovInvXY = cholCovMat^{-1} %*% cbind(obs, covariates)
  cholCovInvXY = Matrix::solve(cholCovMat, obsCov)
#  cholCovInvX = Matrix::solve(cholCovMat, covariates)
  cholCovInvX = cholCovInvXY[,-1]
  # cholCovInvY = cholCovMat^{-1} %*% observations
  #cholCovInvY = Matrix::solve(cholCovMat, observations)
  cholCovInvY = cholCovInvXY[,1]
  
  cholCovInvXcross = Matrix::crossprod(cholCovInvX)
  cholCovInvXcrossInv =       Matrix::solve(cholCovInvXcross)
  detCholCovInvXcross = Matrix::determinant(cholCovInvXcross)$modulus
  
  betaHat = as.vector(
      cholCovInvXcrossInv %*% 
          Matrix::crossprod(cholCovInvX, cholCovInvY)
  ) 
  
  resids = obsCov[,1] - as.vector(obsCov[,-1] %*% betaHat)
  # sigsqhat = resids' %*% Vinv %*% residsx
  #    =   resids' Linv' Linv resids
  cholCovInvResid = Matrix::solve(cholCovMat, resids)
  detCholCovMat = determinant(cholCovMat)$modulus
  
  Nobs = nrow(obsCov)
  Ncov = ncol(obsCov)-1
  Nadj = c(ml=Nobs, reml=Nobs-Ncov)
  totalSsq = as.vector(Matrix::crossprod(cholCovInvResid))
  totalVarHat = totalSsq/Nadj
  
  minusTwoLogLik =  
      Nadj * log(2*pi) +
      2*detCholCovMat  

  minusTwoLogLik = cbind(
      fixedVariance = minusTwoLogLik+totalSsq,
      estimatedVariance = 
           minusTwoLogLik + Nadj + Nadj * log(totalVarHat)
  )
  
  minusTwoLogLik['reml',] = minusTwoLogLik['reml',] + detCholCovInvXcross
  
  varMle = outer(totalVarHat, param[c("variance","nugget")])
  
  names(betaHat) = colnames(cholCovInvXY)[-1]
  cholCovInvXcrossInv=as.matrix(cholCovInvXcrossInv)
  dimnames(cholCovInvXcrossInv)=list(names(betaHat),names(betaHat))
  
  varBetaHat = array(NA, 
      c(dim(cholCovInvXcrossInv),2,2),
      dimnames = c(dimnames(cholCovInvXcrossInv),
          list(
              c('fixedVariance','estimatedVariance'), 
              c('ml','reml')
          )
      )
  )
  
  varBetaHat[,,'fixedVariance','ml'] =
      varBetaHat[,,'fixedVariance','reml'] =
      cholCovInvXcrossInv
  
  varBetaHat[,,'estimatedVariance','ml'] =
      cholCovInvXcrossInv*totalVarHat['ml']
  varBetaHat[,,'estimatedVariance','reml'] =
      cholCovInvXcrossInv*totalVarHat['reml']

  resultR = list(
          minusTwoLogLik = minusTwoLogLik,
          varMle = varMle,
          betaHat = betaHat,
          varBetaHat = varBetaHat,
          totalVarHat = totalVarHat,
          resids=resids,
          detCholCovInvX=detCholCovInvXcross
          )
  } else { # use C
    Nobs = nrow(obsCov)
    Ncov = ncol(obsCov)-1
    Nrep = 1
    
    temp = resultC = .C('maternLogLGivenChol',
     obsCov = as.double(obsCov),
     N= as.integer(c(Nobs,Nrep,Ncov)),
     cholCovMat = as.double(cholCovMat),
     totalSsq = as.double(-9.9),
     betaHat = as.double(rep(-9.9, Ncov)), 
     varBetaHat = as.double(rep(-9.9, Ncov* Ncov)),
     determinants=as.double(c(-9.9,-9.9))
  )
  resultC$detCholCovInvXcrossHalf = resultC$determinants[2]

  Nadj = c(ml=Nobs, reml=Nobs-Ncov)
  totalSsq = resultC$totalSsq
  totalVarHat = totalSsq/Nadj
  
  minusTwoLogLik =  
      Nadj * log(2*pi) +
      2*detCholCovMat
  
  minusTwoLogLik = cbind(
      fixedVariance = minusTwoLogLik+ totalSsq,
      estimatedVariance = 
          minusTwoLogLik + Nadj + Nadj * log(totalVarHat)
  )
  
  minusTwoLogLik['reml',] = minusTwoLogLik['reml',] +
      2*resultC$detCholCovInvXcrossHalf
  
  varMle = outer(totalVarHat, param[c("variance","nugget")])
  
  betaHat = resultC$betaHat
  names(betaHat) = colnames(obsCov)[-1]
  
  varBetaHat = new("dsyMatrix", 
      Dim = as.integer(c(Ncov, Ncov)), 
      uplo="L",
      x=resultC$varBetaHat)
  
  varBetaHat = array(as.matrix(varBetaHat), 
      c(length(betaHat),length(betaHat),2,2),
      dimnames = list(
              names(betaHat), names(betaHat),
              c('fixedVariance','estimatedVariance'), 
              c('ml','reml')
          )
      )

      varBetaHat[,,'estimatedVariance','ml'] =
      varBetaHat[,,'estimatedVariance','ml']*totalVarHat['ml']
  varBetaHat[,,'estimatedVariance','reml'] =
      varBetaHat[,,'estimatedVariance','reml']*totalVarHat['reml']
  
  resultR = list(
      minusTwoLogLik = minusTwoLogLik,
      varMle = varMle,
      betaHat = betaHat,
      varBetaHat = varBetaHat,
      totalVarHat = totalVarHat,
      detCholCovInvX=2*resultC$detCholCovInvXcrossHalf
  )
  

  
  } 
  

  return(resultR)
}