sim.mulmodII<-function(n=NULL,seed=1,N=NULL)
{
d<-2
mixnum<-3
D<-4
M<-matrix(0,mixnum,d)
M[1,]<-c(0,0)   
M[2,]<-c(D,0)  
M[3,]<-c(D/2,D*sqrt(3)/2)   
sig<-matrix(1,mixnum,d)
p<-c(.2,.35,.45)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed)
   return(dendat)
}

if (!is.null(N)){
    #eg<-evalgrid(M,sig,p,N)
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


