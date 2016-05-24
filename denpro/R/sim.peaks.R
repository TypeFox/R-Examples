sim.peaks<-function(n=NULL,seed=1,N=NULL)
{
d<-2
mnum<-5
D<-3.5  
shift<-0.6
M<-matrix(0,mnum,d)
M[1,]<-c(0,0)
M[2,]<-c(D,0)
M[3,]<-c(D/2,D)       #c(D/2,D*sqrt(3)/2)
M[4,]<-c(shift+D/2,D)
M[5,]<-c(-shift+D/2,D)

sig<-matrix(1,mnum,d)
std<-0.3
sig[4,]<-c(std,std)
sig[5,]<-c(std,std)

#p<-c(.25,.35,.45)
p<-c(6/8/3,6/8/3,6/8/3,1/8,1/8)
p<-p/sum(p)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed)
   return(dendat)
}

if (!is.null(N)){
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


