sim.penta4d<-function(n=NULL,seed=1,N=NULL,dist=4)
{
d<-4
moodi<-5
M<-matrix(0,moodi,d)
#dist<-4     # determine the distance between vertices of the pentahedron
M[1,]<-dist*c(1/2, 0,0,0)
M[2,]<-dist*c(-1/2,0,0,0)
M[3,]<-dist*c(0,sqrt(3)/2,0,0)
M[4,]<-dist*c(0,1/(2*sqrt(3)),sqrt(2/3),0)
M[5,]<-dist*c(0,1/(2*sqrt(3)),1/(2*sqrt(6)),sqrt(15/24))
sig<-matrix(1,moodi,d)
p0<-1/moodi
p<-p0*rep(1,moodi)

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


