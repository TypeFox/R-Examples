sim.cross<-function(n=NULL,seed=1,N=NULL,sig1=0.5,sig2=1.5)
{
d<-2
mixnum<-2
M<-matrix(0,mixnum,d)
M[1,]<-c(0,0)      
M[2,]<-c(0,0)      
sig<-matrix(1,mixnum,d)
sig[1,1]<-sig1 
sig[1,2]<-sig2   
sig[2,1]<-sig2   
sig[2,2]<-sig1  
p<-matrix(1,mixnum,1)
p<-p/sum(p)

if (!is.null(n)){
   dendat<-simmix(n,M,sig,p,seed=seed)  
   theta<-pi/4
   rotmat<-matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2)
   dendat<-dendat%*%rotmat
   return(dendat)
}

if (!is.null(N)){
    theta<-pi/4
    eg<-pcf.func("mixt",N,sig=sig,M=M,p=p,theta=theta)
    return(eg)
}

if (is.null(N) && is.null(n)){
   return(list(M=M,sig=sig,p=p,theta=theta))
}

}


