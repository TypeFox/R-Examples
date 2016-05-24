
#Function that computes the BLOM's estimators of the two parameters of Weibull 
BLOMEst<-function(x){
  
  x <- log(x)
  y<- sort(x)
  n<- length(y)
  I<- 1:1:n
  Theta<- c(0,-(1-I/(n+1))*log(1-I/(n+1)),0)
  C1<- Theta[1:(n+1)] -Theta[2:(n+2)] 
  Esp0 <- GoFNS(1,n,n)
  Esp<- c(0,Esp0,0)
  C2<- Theta[1:(n+1)]*Esp[1:(n+1)]-Theta[2:(n+2)]*Esp[2:(n+2)]
  mat <- matrix(c(sum(C1^2),sum(C1*C2),sum(C1*C2),sum(C2^2)),nrow=2,ncol=2)
  det <-  determinant(mat,logarithm=FALSE)$modulus
  I11 = 1/det*mat[2,2]
  I21 = -1/det*mat[2,1]
  I12 = -1/det*mat[1,2]
  I22 = 1/det*mat[1,1] 
  g1=Theta[2:(n+1)]*{I11*(C1[2:(n+1)]-C1[1:n])+I12*(C2[2:(n+1)]-C2[1:n])}
  g2=Theta[2:(n+1)] *{I21*(C1[2:(n+1)]-C1[1:n])+I22*(C2[2:(n+1)]-C2[1:n])}  
  ksi = sum(g1*y)
  t=sum(g2*y)
  #Compute the pseudo-observation y_1, .., y_n
  y = (y-ksi)/t
  return(BLOMEst<-list(eta=exp(ksi),beta=1/t,y=y))
}
