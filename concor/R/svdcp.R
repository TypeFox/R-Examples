"svdcp" <-
function(x,H,r) { 
# Initialisations 
q<-dim(x)[2]
if (sum(H) != q) print("YOUR H IS NOT SUITABLE")

k<-length(H)
s2<-matrix(0,k,r)
u<-matrix(0,dim(x)[1],r);
v<-matrix(0,q,r);
kx<-cumsum(c(0,H));

# Calculus
 for (i in 1:r) {
  s<-svd(x)
  u[,i]<-s$u[,1]
    for (j in 1:k) {
ax <- (kx[j]+1):kx[j+1]
norm2 <- t(s$v[ax,1]) %*% s$v[ax,1]
      s2[j,i]<-norm2 * s$d[1]^2 
      if (s2[j,i] > 1e-8) {
v[ax,i]<-s$v[ax,1]/sqrt(norm2)
x[,ax] <- x[,ax]-x[,ax]%*%(v[ax,i]%*%t(v[ax,i]))
       }

    }
 }

list(u=u,v=v,s2=s2)
}

