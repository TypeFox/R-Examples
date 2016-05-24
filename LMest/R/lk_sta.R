lk_sta <- function(tau,u,V,G2,outl=TRUE){

  k = length(u)
  PI = exp(G2%*%tau)
  PI = t(matrix(PI,k,k))
  PI = diag(1/rowSums(PI))%*%PI
  #la = colMeans(PI%^%10000) 
  PI1 = PI
  for (i in 1:10000) PI1 = PI1%*%PI
  la = colMeans(PI1)
  la = la/sum(la)
  lk = u%*%log(la)+sum(V*log(PI))
  flk = -lk
  if (outl) out = list(flk=flk,la=la,PI=PI)	else flk	
}


