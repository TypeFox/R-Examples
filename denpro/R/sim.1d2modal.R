sim.1d2modal<-function(n=NULL,seed=1,N=NULL,distr=FALSE)
{
d<-1
M<-c(0,2,4)
mixnum<-length(M)
sig<-matrix(1,mixnum,d)
sig[1]<-0.3
p<-matrix(1,mixnum,1)
p[2]<-2
p<-p/sum(p)

if (!is.null(n)){
   dendat<-simmix(n=n,M,sig,p,seed=seed,d=1)
   return(dendat)
}

if (!is.null(N)){
    xala<--2
    xyla<-7
    support<-c(xala,xyla)
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p,support=support,distr=distr)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p))
}

}


