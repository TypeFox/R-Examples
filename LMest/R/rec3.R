rec3 <- function(Q,PI,Pio,pim){

# preliminaries
  out = dim(Pio)
  n = out[1]; k = out[2]; TT = out[3] 
# for t=T
  if(k==1){
    U = array(1,c(k,n,TT)); V = matrix(n*(TT-1),k,k)
  }else{
    U = array(0,c(k,n,TT)); V = matrix(0,k,k)
    U[,,TT] = t(Q[,,TT]/pim)
    Q1 = t(Q[,,TT-1]/pim)
    P1 = Pio[,,TT]
    R = matrix(1,k,n);
    for(i in 1:n) V = V+Q1[,i]%o%P1[i,]
# for other t
    for(t in seq(TT-1,2,-1)){
      R = PI%*%(t(Pio[,,t+1])*R)
      P1 = Pio[,,t]*t(R)
      U[,,t] = Q1*R
      Q1 = t(Q[,,t-1]/pim)
      for(i in 1:n) V = V+Q1[,i]%o%P1[i,]
    }
    V = PI*V
    R = PI%*%(t(Pio[,,2])*R)
    U[,,1] = Q1*R
  }
  out = list(U=U,V=V)
  out
}