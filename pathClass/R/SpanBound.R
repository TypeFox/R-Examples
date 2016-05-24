spanbound <- function(fit, xtrn, ytrn){
  
  svindex = unlist(alphaindex(fit))	
  alpha = unlist(coef(fit))	
  pos = which(alpha > 0)
  neg = which(alpha < 0)
  alpha = abs(alpha)	
  if(class(xtrn) != "kernelMatrix"){
    yhat = predict(fit, xtrn, type="decision")
    if(param(fit)$C != Inf)
      error("span bound is only for L2-SVM!")
    K = kernelMatrix(kernelf(fit), xtrn)
  }
  else{
    K = xtrn
    yhat = K[,svindex]%*%as.matrix(unlist(coef(fit))) - b(fit)
  }
  output = ytrn*yhat	
  Cpos = Inf
  Cneg = Inf
  eps = 1e-5	
  boundpos = (alpha[pos] >= Cpos*(1-eps))
  boundneg = (alpha[neg] >= Cneg*(1-eps))
  sv1pos = svindex[pos[!boundpos]]
  sv2pos = svindex[pos[boundpos]]
  sv1neg = svindex[neg[!boundneg]]
  sv2neg = svindex[neg[boundneg]]	
  sv1 = sort(c(sv1pos, sv1neg))
  sv2 = sort(c(sv2pos, sv2neg))
  n = ncol(K)
  span = double(n)
  alpha1 = double(n)
  alpha1[svindex] = alpha 		
  if(length(sv1) > 0){ # in-bound SVs 								
    ell = length(sv1)	
    invK = chol2inv(chol(K[sv1,sv1,drop=FALSE]))
    T = -1/sum(invK)
    T2 = invK%*%as.matrix(rep(1,ell))*T
    T3 = t(as.matrix(rep(1,ell)))%*%invK
    invKSV = rbind(cbind(invK + T2%*%T3, -T2), cbind(-T*T3, T))
    tmp = diag(as.matrix(invKSV)) + 1e-10
    span[sv1] = 1./tmp[1:ell]
  }	
  else
    warning("No in-bound SVs!")			
  if(length(sv2) > 0){	# bound SVs	
    span[sv2] = diag(as.matrix(K[sv2,sv2,drop=FALSE]))
    if(length(sv1) > 0){
      V = rbind(K[sv1,sv2,drop=FALSE], rep(1,length(sv2)))			
      span[sv2] = span[sv2] - diag(t(V)%*%invKSV%*%V)			
    }				 
  }	
  loo = mean((output - alpha1*span <= 0)*1)				
  ## cat("Span bound =", loo,"\n\n")		
  loo
}
