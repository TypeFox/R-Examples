`hbrr.integrate3.rhoeq1` <-
function(MU,V,A=c(1,1,1),...){
    if (V[1,1]!=V[2,2] | V[2,2]!=V[3,3] | MU[1]!=MU[2] | MU[2]!=MU[3]){
        stop("when rho=1, then must have V[1,1]=V[2,2]=V[3,3], and MU[1]=MU[2]=MU[3]") } 
    hb.x<-function(x){
            out<-(1/(1+10^(A[1]*x)))*(1/(1+10^(A[2]*x)))**(1/(1+10^(A[3]*x)))*dnorm(x,mean=MU[1],sd=sqrt(V[1,1]))
            return(out)
        }
   out<-integrate(hb.x,-Inf,Inf,...)$value
   return(out)
}

