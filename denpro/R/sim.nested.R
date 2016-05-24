sim.nested<-function(n=NULL,seed=1,N=NULL)
{
d<-2
con<-3^(3/2)/2
shift<-0.8
std<-0.3
mnum<-5
M<-matrix(0,mnum,d)
M[1,]<-c(0,con)
M[2,]<-c(3,-con)
M[3,]<-c(-3,-con) 
M[4,]<-c(shift,con)
M[5,]<-c(-shift,con)

sig<-matrix(1,mnum,d)
sig[4,]<-c(std,std)
sig[5,]<-c(std,std)
p<-c(5/7/3,5/7/3,5/7/3,1/7,1/7)
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


