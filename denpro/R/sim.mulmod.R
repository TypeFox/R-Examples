sim.mulmod<-function(n=NULL,seed=1,N=NULL)
{
d<-2
mnum<-3
D<-3.5  #3 
shift<-0.2  
M<-matrix(0,mnum,d)
M[1,]<-c(0,0)
M[2,]<-c(D,shift)   #c(D,0)
M[3,]<-c(D/2,D)    #   #c(D/2,D*sqrt(3)/2)
sig<-matrix(1,mnum,d)
p<-c(.25,.35,.45)
p<-p/sum(p)

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


