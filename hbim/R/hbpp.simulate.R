`hbpp.simulate` <-
function(MU,V,A,RP=.1,NSIM=10^5){
    nmu<-length(MU)
    if (nmu==1){
         x<-rnorm(NSIM,MU,sd=sqrt(V))
         protected<-rep(1,NSIM)
         protected[(1/(1+10^(A*x)))>RP]<-0
         out<-mean(protected)
    }
    else{
        fprot<-function(d){
            RR<-prod(1/(1+10^(A*d)))
            if (RR>RP){ protected<-0 }
            else{ protected<-1 }
            return(protected)  
        }
        x<-rmvnorm(NSIM,mean=MU,sigma=V)
        out<-mean(apply(x,1,fprot))
    }
    out
}

