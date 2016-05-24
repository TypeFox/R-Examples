pMaxCorrNor<-function(x,k,rho){
    inner.int<-function(s){
      pnorm((s*sqrt(rho)+x)/sqrt(1-rho))^k*dnorm(s)
    }
    return(1-integrate(inner.int,-Inf,Inf)$value)
}