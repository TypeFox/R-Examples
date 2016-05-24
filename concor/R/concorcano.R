"concorcano" <-
function(x,y,py,r) {
# INITIALISATIONS 

n<-dim(x)[1]
q<-dim(y)[2]
if (sum(py) != q ) stop("py IS NOT SUITABLE")

s<-svd(x);rk<-sum(s$d > max(dim(x))*s$d[1]*1e-8)
P<-matrix(s$u[,1:rk]*sqrt(n),ncol=rk)

ky<-length(py)
ry<-matrix(0,1,ky)

Py<-NULL
cuy=c(0,cumsum(py))
for (j in 1:ky) {
 s<-svd(y[,(cuy[j]+1):cuy[j+1]])
 ry[j]<-sum(s$d > max(c(n,py[j]))*s$d[1]*1e-8)
 Py<-cbind(Py,s$u[,1:ry[j]])
}

if (r > min(c(min(ry),rk,n))) stop("r IS TOO HIGH")
Py<-matrix(Py,ncol=sum(ry))*sqrt(n)
s<-concor(P,Py,ry,r);
cumuly<-cumsum(c(0,ry))
cy<-matrix(0,ky*n,r)
   for (k in 1:ky) {
ay<-(cumuly[k]+1):cumuly[k+1]
cy[(n*(k-1)+1):(n*k),]<-matrix(Py[,ay],ncol=ry[k])%*%s$v[ay,]
   }
list(cx=P%*%s$u,cy=cy,rho2=s$cov2)
}

