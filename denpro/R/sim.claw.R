sim.claw<-function(n=NULL,seed=1,N=NULL)
{
d<-1

M<-c(0,-1,-0.5,0,0.5,1)
sig<-c(1,0.1,0.1,0.1,0.1,0.1)
p<-c(0.5,0.1,0.1,0.1,0.1,0.1)
mixnum<-length(M)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed,d=1)
   return(dendat)
}

if (!is.null(N)){
    support<-c(-3,3)
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p,support=support)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


