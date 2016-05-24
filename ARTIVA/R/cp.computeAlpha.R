cp.computeAlpha <-
function(birth,lNew,kminus,Ekl,Estar,Ekr,yL,PxL,yR,PxR,y2,Px2,D,delta2, q, smax, v0, gamma0){
  # birth = 1 for birth, -1 for death.
  # lNew = number of edges in the new phase
  # kminus = minimal number of changepoints between the 2 compared models (=s for birth, s-1 for death) -> INUTILE !!!
  # Ekl = 
  # Estar =
  # yL,yR, y2 : response data  (left, right, both)
  # PxL, PxR, Px2 : projection matrix (left, right, both)
  # D : hyperparms for the number of edges in each phase. (Number of edges s ~ truncated Poisson P(D). )
  # delta2 : hyperparms for empirical covariance (can be seen as the expected  signal-to-noise ratio)  ~IG(alphad2,betad2)
  #q
  #smax
  #v0
  #gamma0
  
  # last modified by Sophie 16/09/09 
  logR=  +(v0/2)*log(gamma0/2)-lgamma(v0/2)  -log( (sqrt(delta2+1))^(lNew+1) ) +lgamma(((Estar-Ekl)+v0)/2)+lgamma(((Ekr-Estar)+v0)/2)-lgamma(((Ekr-Ekl)+v0)/2) +(((Estar-Ekl)+v0)/2)*log(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))  +(((Ekr-Estar))/2)*log (((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))  -(v0/2)*log((gamma0+t(yR) %*% PxR %*% yR)/2)

  logR=birth*logR

  if(logR>0){
    res=1
  } else {
    res=exp(logR)
  }
  
  return(res)
}
