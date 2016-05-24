"concor" <-
function(x,y,py,r) {
# INITIALISATIONS 
n<-dim(x)[1]
p<-dim(x)[2]
q<-dim(y)[2]

if (sum(py) != q ) stop("py IS NOT SUITABLE")
if (r > min(c(min(py),q,n))) stop("r IS TOO HIGH")

ky<-length(py)
cri<-matrix(0,ky,r)
cumuly=cumsum(c(0,py))
u<-matrix(0,p,r)
V<-matrix(0,q,r)
v<-V

for (i in 1:r) {
  s<-svd(t(x)%*%y)
  u[,i]<-s$u[,1]
  V[,i]<-s$v[,1]
  c1=s$d[1]^2
   for (k in 1:ky) {
ay<-(cumuly[k]+1):cumuly[k+1]
ny<-t(V[ay,i])%*%V[ay,i]
cri[k,i]<-ny*c1 
        if (ny > 1e-8) {
 v[ay,i]<-V[ay,i]/sqrt(ny)
         y[,ay]<-y[,ay]-y[,ay]%*%(v[ay,i]%*%t(v[ay,i]))
        }
   }
}
list(u=u,v=v,V=V,cov2=cri/n^2)
}

