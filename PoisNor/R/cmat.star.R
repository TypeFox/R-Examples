cmat.star <-
function(no.pois, no.norm, corMat, lamvec){
  
  nPois = length(lamvec)
  nNorm = ncol(corMat) - nPois
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  Cor_NNforPP<-function(lambda1, lambda2, r){
    samples=100000
    u = runif(samples, 0, 1)
    lambda=c(lambda1,lambda2)
    maxcor=cor(qpois(u, lambda1), qpois(u, lambda2))
    mincor=cor(qpois(u, lambda1), qpois(1-u, lambda2))
    a=-maxcor*mincor/(maxcor+mincor)
    b=log((maxcor+a)/a, exp(1))
    c=-a
    corrected=log((r+a)/a, exp(1))/b
    corrected=ifelse ((corrected>1 | corrected<(-1)),NA, corrected)
    return(corrected)
  }
  
  #-------------------------------------------------------------
  Cor_NNforPN<-function(PN.cor, lam){
    X=rnorm(100000,0,1)
    Y=rnorm(100000,0,1)
    
    U = pnorm(X)
    
    Xpois = qpois(U,lam)
    
    c = cor(Xpois[order(Xpois)],Y[order(Y)])/cor(X[order(X)],Y[order(Y)]) 
    r = PN.cor/c 
    return(r)
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  

    if (Validate.correlation (no.pois, no.norm, corMat, lamvec) ){

    corMat.star = diag(nrow(corMat) )
    
      for (i in 1:nrow(corMat)){
        for (j in 1:nPois){
          
          if (i!=j & i <=nPois ){ 
            corMat.star [i,j] = Cor_NNforPP(lamvec[i], lamvec[j], corMat[i,j] )
          }
          
          if (i > nPois) {
            corMat.star [i,j] = corMat.star [j,i] = Cor_NNforPN(corMat[i,j], lamvec[j] )
          }
          cat(".")
        }
      }
      if (no.norm>0){
        corMat.star  [(nPois +1): nrow(corMat),  (nPois +1): ncol(corMat)] = corMat[(nPois +1): nrow(corMat),  (nPois +1): ncol(corMat)]
      }
    }
  cat("\n")
  
  if(!is.positive.definite(corMat.star)){
    warning( "Intermediate correlation matrix is not positive definite. Nearest positive definite matrix is used!")
    corMat.star=as.matrix(nearPD(corMat.star, corr = TRUE, keepDiag = TRUE)$mat)
  }
  corMat.star = ( corMat.star+t(corMat.star) )/2
  
  return (corMat.star)
    
}
