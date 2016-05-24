`hbrr.integrate3` <-
function(MU,V,A=c(1,1,1),...){
    hb.z<-function(z){
        hb.yZ<-function(y,Z){
            ny<-length(y)
            outy<-rep(NA,ny)
            hb.xYZ<-function(x,Y,Z){
                q<-matrix(c(x,rep(Y,length(x)),rep(Z,length(x))),ncol=3)
                outx<-(1/(1+10^(A[1]*x)))*
                      (1/(1+10^(A[2]*Y)))*
                      (1/(1+10^(A[3]*Z)))*
                      dmvnorm(q,mean=MU,sigma=V)
                return(outx)
            }
            for (j in 1:ny){
                ### integrate over x keeping y and z fixed
                outy[j]<-integrate(hb.xYZ,-Inf,Inf,Y=y[j],Z=Z,...)$value
            }
            outy
       }
       nz<-length(z)
       outz<-rep(NA,nz)
       for (k in 1:nz){
           outz[k]<-integrate(hb.yZ,-Inf,Inf,Z=z[k],...)$value
       }
       return(outz)
    }
    out<-integrate(hb.z,-Inf,Inf,...)$value
    return(out)
}

