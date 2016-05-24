sim.fox<-function(n=NULL,seed=1,N=NULL)
{
d<-2
mixnum<-14
D<-1.8
M<-matrix(0,mixnum,d)
M[1,]<-c(0,0)      #c(0,0)

M[2,]<-c(D,0)      #c(D1,0)
M[3,]<-c(2*D,0)

M[4,]<-c(0,D)
M[5,]<-c(0,2*D)
M[6,]<-c(0,3*D)

M[7,]<-c(0,-D)
M[8,]<-c(0,-2*D)
M[9,]<-c(0,-3*D)

M[10,]<-c(1.5,3.9*D)
M[11,]<-c(-1.5,3.7*D)
M[12,]<-c(-1.5,4.2*D)
M[13,]<-c(-1.5,4.5*D)
M[14,]<-c(-1.5,4.7*D)

sig<-matrix(1,mixnum,d)
sig[10,1]<-0.7
sig[11,1]<-0.7
sig[12,1]<-0.7
sig[13,1]<-0.7
sig[14,1]<-0.7
p<-matrix(1,mixnum,1)
p[6]<-0.6
p[10]<-0.3
p[11]<-0.25
p[12]<-0.1
p[13]<-0.05
p[14]<-0.05
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


