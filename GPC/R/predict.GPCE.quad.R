predict.GPCE.quad <- function(object,n,x,...) {
    
  p <- object$Args$p
  d <- object$Args$InputDim
  m = getM(d,p)
  Index = indexCardinal(d,max(p))  
  PCECoeff <- object$PCEcoeff
  
  if (is.null(object$Args$ParamDistrib$alpha)){ alpha <- NULL    
  } else { alpha <- object$Args$ParamDistrib$alpha }
  if (is.null(object$Args$ParamDistrib$beta)){ beta <- NULL    
  } else { beta <- object$Args$ParamDistrib$beta }

  y <- rep(0,n)
  for (nn in 1:n){
    for (mm in 1:m){
      tmp = 1
      for (dd in 1:d){
        tmp = tmp*polyNorm(Index[dd,mm],x[nn,dd],object$Args$ExpPoly[dd],alpha[dd],beta[dd])
      }      
      y[nn] = y[nn] + PCECoeff[mm]*tmp
    }
  }
  return(y)
}
  #   # Poly <- hermite.he.polynomials(p, normalized=FALSE)
  #   y <- rep(0,n)
  #   for (nn in seq(1,n)){
  #     for (mm in seq(1,getM(d,p))){
  #       tmp <- 1;
  #       for (dd in seq(1,d)){
  #         polyNorm(Index[nn,mm],NDXi[nn,zz],PolyName[nn],alpha,beta)
  #         tmp = tmp * unlist(polynomial.values(Poly[index[dd,mm]+1],x[dd,nn]))
  #       }     
  #       y[nn] = y[nn] + PCECoeff[mm]*tmp
  #     }
  #   }
  #   return(y)
  #   A=object$TruncSet[Selection,]
  #   Design=GetDesign4Predict(x,object$Args$InputDistrib,object$Args$ParamDistrib,object$Args$PCSpace)
  #   DM=DataMatrix(A,Design$PCE,object$Args$InputDistrib,object$Args$pmaxi,object$Args$PCSpace)
  #   MetaOut=DM%*%object$CoeffPCE[Selection,]  
  #   return(MetaOut)