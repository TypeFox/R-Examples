speclagreg.K.threshold = function(X,Y,freq){  
  K = freq
  n = dim(X)[1]
  th = 1/sqrt(n)
  
  SXX = spectral.density(X,freq=freq)
  for (i in 1:length(K)){
    E = eigen(SXX$operators[i,,])
    K[i] = max(1,sum(abs(E$values) > th))
  }
  
  return(K)    
}  
