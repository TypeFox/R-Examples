`hbrr.integrate1` <-
function(MU,V,A=1,...){
    hb.x<-function(x){
            out<-(1/(1+10^(A*x)))*dnorm(x,mean=MU,sd=sqrt(V))
            return(out)
        }
   out<-integrate(hb.x,-Inf,Inf,...)$value
   return(out)
}

