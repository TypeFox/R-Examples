# Method 'saddlepoint' from package 'survey' v.3.30-3 (c) 2014 Thomas Lumley

pchisqsum<-function(x,df,a,lower.tail=TRUE, method="saddlepoint"){

  satterthwaite<-function(a,df){
    if(any(df>1)){
      a<-rep(a,df)
    }
    tr<-mean(a)
    tr2<-mean(a^2)/(tr^2)
    
    list(scale=tr*tr2, df=length(a)/tr2)
  }
  
  sat<-satterthwaite(a,df)
  guess<-pchisq(x/sat$scale,sat$df,lower.tail=lower.tail)

  for(i in seq(length=length(x))){
    lambda<-rep(a,df)
    sad<-sapply(x,saddle,lambda=lambda)
    if (lower.tail) sad<-1-sad
    guess<-ifelse(is.na(sad),guess,sad)
  }

  return(guess)

}

saddle<-function(x,lambda){
  d<-max(lambda)
  lambda<-lambda/d
  x<-x/d
  k0<-function(zeta) -sum(log(1-2*zeta*lambda))/2
  kprime0<-function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
  kpprime0<-function(zeta) 2*sum(lambda^2/(1-2*zeta*lambda)^2)
  n<-length(lambda)
  if (any(lambda < 0)) {
    lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
  } else if (x>sum(lambda)){
    lmin <- -0.01
  } else {
    lmin<- -length(lambda)/(2*x)
  }
  lmax<-min(1/(2*lambda[lambda>0]))*0.99999

  hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, 
                       lower = lmin, upper = lmax, tol = 1e-08)$root
  
  w<-sign(hatzeta)*sqrt(2*(hatzeta*x-k0(hatzeta)))
  v<-hatzeta*sqrt(kpprime0(hatzeta))
  if (abs(hatzeta)<1e-4)
    NA
  else
    pnorm(w+log(v/w)/w, lower.tail=FALSE)
}

