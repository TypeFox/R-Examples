#Function that computes the maximum likelihood estimators of the two parameters of Weibull 
MLEst<-function(x){
  if(sum(x<0)){stop(paste("Data x is not a positive sample"))}
  n=length(x)
  #Apply the log. transformation => If: X -> Weibull, then -log(X) -> EV
  lv = -log(x)
  y = sort(lv)
  #Compute the MLEs of the EV distribution
  f1 <- function(toto,vect){
    if(toto!=0){
      f1= sum(vect)/length(vect)
      f1 = f1 - sum(vect*exp(-vect*toto))/sum(exp(-vect*toto)) - 1/toto
    }else{ f1=100}
    f1=abs(f1)
  }  
  theta <- optimize(f1,c(0.0001,50),maximum=FALSE,vect=y,tol=1/10^5)
  t <- theta$minimum
  aux <- sum(exp(-y*t)/length(x))
  ksi <- -(1/t)*log(aux)
  #Compute the pseudo-observation y_1, .., y_n
  y <- -(y-ksi)*t
  return(MLEst<-list(eta=exp(-ksi),beta=t,y=y))
}

#Function that computes the maximum likelihood estimators of the two parameters of Weibull in the case of right simple censoring type II
 MLEst_c<-function(x,r){
  if(sum(x<0)){stop(paste("Data x is not a positive sample"))}
  x = sort(x)
  l = length(x)
  n = l+r

 #Compute beta
  f1 <- function(c,vect){
  t = length(vect)
  if(c!=0){
   f1= 1/c+1/(n-r)*sum(log(vect))
   f1 = f1 - sum(log(vect)*vect^c+r*vect[t]^c*log(vect[t]))/sum(vect^c+r*vect[t]^c)
  }else{ f1=100}
   f1=abs(f1)
  }
  op <- optimize(f1,c(0.0001,50),maximum=FALSE,vect=x,tol=1/10^5)
  b <- op$minimum
  eta <- ((sum(x^b)+r*x[l]^b)/(n-r))^(1/b)
  y <- (log(x)-log(eta))*b
  return(MLEst<-list(eta=eta,beta=b,y=sort(y)))

}
#######################################
