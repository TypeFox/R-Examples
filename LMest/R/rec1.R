rec1 <- function(Pio,las,PI){

  out = dim(Pio)  
  n = out[1]; k = out[2]; TT = out[3] 
  La = rep(1,n)%o%las
  Q = array(0,c(n,k,TT))
  Q[,,1] = Pio[,,1]*La
  for(t in 2:TT) Q[,,t] = Pio[,,t]*(Q[,,t-1]%*%PI)
  Q
 
}