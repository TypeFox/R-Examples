`hbrr.simulate` <-
function(MU,V,A,NSIM=10^4){
    nmu<-length(MU)
    if (nmu==1){
         x<-rnorm(NSIM,MU,sd=sqrt(V))
         out<-mean(1/(1+10^(A*x)))
    }
    else{
        hill.bliss<-function(d){
            prod(1/(1+10^(A*d)))
        }
        x<-rmvnorm(NSIM,mean=MU,sigma=V)
        out<-mean(apply(x,1,hill.bliss))
    }
    out
}

