
pchisqsum<-function(x,df,a,lower.tail=TRUE,
                    method=c("satterthwaite","integration","saddlepoint")){
  
  satterthwaite<-function(a,df){
    if(any(df>1)){
      a<-rep(a,df)
    }
    tr<-mean(a)
    tr2<-mean(a^2)/(tr^2)
    
    list(scale=tr*tr2, df=length(a)/tr2)
  }
  
  ## chisqphi<-function(t, df=1,a=1){
  ##     (1-(0+2i)*a*t)^(-df/2)
  ##   }
  
  ## make.integrand<-function(x,DF,A){ 
  ##   m<-length(DF)
    
  ##   function(t){
  ##     n<-length(t)
  ##     tmp<-matrix(chisqphi(rep(t,each=m),rep(DF,n),rep(A,n) ),ncol=n)
  ##     phi<-apply(tmp,2,prod)
  ##     rval<-Im(phi*exp(-(0+1i)*t*x)/(pi*t))
  ##     rval[t==0]<-x/(pi)
  ##     rval
  ##   }
    
  ## }
  
  
  method<-match.arg(method)
  sat<-satterthwaite(a,df)
  guess<-pchisq(x/sat$scale,sat$df,lower.tail=lower.tail)
  
  if (method=="satterthwaite")
    return(guess)
  
  method<-match.arg(method)
  if (method=="integration" && !suppressWarnings(require("CompQuadForm"))){
      warning("Package 'CompQuadForm' not found, using saddlepoint approximation")
      method<-"saddlepoint"
    }
  

  abstol<-guess/1000
  abstol<-pmax(1e-9, abstol)
  reltol<-rep(1/1000,length(abstol))
 
  if (method=="integration"){
    require("CompQuadForm")
    if (any(a<=0)){
      for(i in seq(length=length(x))){
        f<-davies(x[i],a,df,acc=1e-7)
        if (f$ifault>0) warning("Probable loss of accuracy ")
        guess[i]<-f$Qq
      }
      if(any(guess<1e-6)) warning("Probable loss of accuracy ")
    } else{
      for(i in seq(length=length(x))){
        guess[i]<-farebrother(x[i],a,df)$res
      }
      if(any(guess<1e-9)) warning("Probable loss of accuracy ")
    }
    if (lower.tail)
      guess<-1-guess
    ## for(i in seq(length=length(x))){
    ##   if (guess[i]< 1e-5 || guess[i]> 1-1e-5)
    ##     next ## don't even try.
    ##   ff<-make.integrand(x[i],df,a)
    ##   rval<-integrate(ff,0,Inf,subdivisions=10000,
    ##                       abs.tol=abstol[i],rel.tol=reltol[i],stop.on.error=FALSE)
    ##   if (inherits(rval, "try-error") || rval$message!="OK"){
    ##     warning("integration failed for x=",x[i],", using Satterthwaite approximation")
    ##   }else
    ##   guess[i]<- if (lower.tail) 1/2-rval$value else 1/2+rval$value
    ## }
    return(guess)
  } else if (method=="saddlepoint"){
    for(i in seq(length=length(x))){
      lambda<-rep(a,df)
      sad<-sapply(x,saddle,lambda=lambda)
      if (lower.tail) sad<-1-sad
      guess<-ifelse(is.na(sad),guess,sad)
    }
    return(guess)
  }
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

