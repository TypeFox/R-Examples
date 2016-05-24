"concorgmcano" <-
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
cux<-c(0,cumsum(rx))
Px<-matrix(Px,nrow=n)
ky<-length(py)
ry<-matrix(0,1,ky)
Py<-NULL
cuy=c(0,cumsum(py))
 for (j in 1:ky) {
  s<-svd(y[,(cuy[j]+1):cuy[j+1]])
  ry[j]<-sum(s$d > max(c(n,py[j]))*s$d[1]*1e-8)
  Py<-cbind(Py,s$u[,1:ry[j]]*sqrt(n))
 }
if (r > min(c(min(ry),min(rx),n))) stop("r IS TOO HIGH")
cuy<-c(0,cumsum(ry))
Py<-matrix(Py,nrow=n)
s<-concorgm(Px,rx,Py,ry,r)
cy<-matrix(0,n*ky,r)
cx<-matrix(0,n*kx,r)

for  (j in 1:kx) {
cx[((j-1)*n+1):(j*n),]<-matrix(Px[,(cux[j]+1):cux[j+1]],nrow=n)%*%s$u[(cux[j]+1):cux[j+1],]
}
for  (j in 1:ky) cy[((j-1)*n+1):(j*n),]<-matrix(Py[,(cuy[j]+1):cuy[j+1]],nrow=n)%*%s$v[(cuy[j]+1):cuy[j+1],]

list(cx=cx,cy=cy,rho2=s$cov2)
}

