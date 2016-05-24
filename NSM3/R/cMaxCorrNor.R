cMaxCorrNor<-function(alpha,k,rho){
  r<-function(x,n){
    inner.int<-function(s){
      pnorm((s*sqrt(rho)+x)/sqrt(1-rho))^n*dnorm(s)
    }
    return(integrate(inner.int,-Inf,Inf)$value)
  }
  pts<-seq(3.5,0,by=-0.01)
  for(i in 1:length(pts)){
    if(r(pts[i],k)<=(1-alpha)){
      return(pts[i])
    }
  }
  
}
