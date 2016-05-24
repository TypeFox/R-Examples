`hbrr.integrate2.rhoeq1` <-
function(MU,V,A=c(1,1),...){
    if (V[1,1]!=V[2,2] | MU[1]!=MU[2]){
        stop("when rho=1, then V[1,1] must equal V[2,2], and MU[1] must equal MU[2]") } 
    hb.x<-function(x){
            out<-(1/(1+10^(A[1]*x)))*(1/(1+10^(A[2]*x)))*dnorm(x,mean=MU[1],sd=sqrt(V[1,1]))
            return(out)
        }
   out<-integrate(hb.x,-Inf,Inf,...)$value
   return(out)
}

