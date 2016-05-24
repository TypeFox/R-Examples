"concorsreg" <-
function(x,px,y,py,r) { 
if (sum(px) != dim(x)[2] | sum(py) != dim(y)[2] ) stop("px or py IS NOT SUITABLE")
n<-dim(x)[1]
kx<-length(px)
rx<-matrix(0,1,kx)
Px<-NULL
cux=c(0,cumsum(px))
 for (j in 1:kx) {
  s<-svd(x[,(cux[j]+1):cux[j+1]])
  rx[j]<-sum(s$d > max(c(n,px[j]))*s$d[1]*1e-8)
  Px<-cbind(Px,s$u[,1:rx[j]]*sqrt(n))
 }
if (r > min(c(min(py),min(rx),n))) stop("r IS TOO HIGH")

cux<-c(0,cumsum(rx))
Px<-matrix(Px,ncol=cux[kx+1])
s<-concors(Px,rx,y,py,r)
cx<-matrix(0,n*kx,r)
for  (j in 1:kx) cx[((j-1)*n+1):(j*n),]<-matrix(Px[,(cux[j]+1):cux[j+1]],nrow=n)%*%s$u[(cux[j]+1):cux[j+1],]
list(cx=cx,v=s$v,varexp=s$cov2)
}

