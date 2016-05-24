sim.tetra3d<-function(n=NULL,seed=1,N=NULL)
{
dist<-3    # determine the distance between vertices of the tetrahedron

d<-3
moodi<-4
M<-matrix(0,moodi,d)
height<-sqrt(3)/2           # sqrt(3)/2 = 0.8660254
len<-1/(2*sqrt(3))          # 1/(2*sqrt(3)) = 0.2886751
kor<-sqrt(2/3)              # sqrt(2/3) = 0.8164966
M[1,]<-dist*c(1/2,0,0)      # ( 1.5, 0.0, 0.0)
M[2,]<-dist*c(-1/2,0,0)     # (-1.5, 0.0, 0.0)
M[3,]<-dist*c(0,height,0)   # ( 0.0, 2.6, 0.0)
M[4,]<-dist*c(0,len,kor)    # ( 0.0, 0.9, 2.4)
sig<-matrix(1,moodi,d)
p0<-1/moodi
p<-p0*rep(1,moodi)

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


