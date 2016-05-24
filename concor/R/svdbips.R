"svdbips" <-
function(x,K,H,r) {
# INITIALISATIONS 
p<-dim(x)[1]
q<-dim(x)[2]
if (sum(H) != q | sum(K) != p) stop("K or H IS NOT SUITABLE")
if (r > min(c(K,H))) stop("r IS NOT SUITABLE")
M<-length(K)
N<-length(H)
u<-matrix(0,p,r)
v<-matrix(0,q,r)
ck<-cumsum(c(0,K))
ch<-cumsum(c(0,H))
s2<-array(0,c(M,N,r))

 #PROPOSED INITIALISATION OF THE ALGORITHM with u and v
    for (i in 1:M) {
      ak<-(ck[i]+1):ck[i+1]
      s<-svd(matrix(x[ak,],nrow=length(ak)))
      u[ak,]<-s$u[,1:r]
    }

    for (j in 1:N) {
      ah<-(ch[j]+1):ch[j+1]
      s<-svd(matrix(x[,ah],ncol=length(ah)))
      v[ah,]<-s$v[,1:r]
    }
cc<-2;cc1<-0

#ALGORITHM
while (abs(cc-cc1) > 1e-8) {
  #aa and bb are converging to the optimized criterion
 aa=0;bb=0;
 cc1=cc;
 A<-matrix(0,p,r)
 B<-matrix(0,r,q)

for (i in 1:M) {
ak<-(ck[i]+1):ck[i+1]
  for (j in 1:N) {
          ah<-(ch[j]+1):ch[j+1]
  d<-diag(t(u[ak,])%*%x[ak,ah]%*%v[ah,]);l<-length(d)
  A[ak,]<-A[ak,]+matrix(x[ak,ah],nrow=K[i])%*%v[ah,]%*%diag(d,nrow=l)
  }
s<-svd(A[ak,]);u[ak,]<-s$u[,1:r]%*%t(s$v)
        aa<-aa+sum(s$d)
}

  for (j in 1:N) {
ah<-(ch[j]+1):ch[j+1]
  for (i in 1:M) {
          ak<-(ck[i]+1):ck[i+1]
          d<-diag(t(u[ak,])%*%x[ak,ah]%*%v[ah,]);l<-length(d)
  B[,ah]<-B[,ah]+diag(d,nrow=l)%*%t(u[ak,])%*%x[ak,ah]
  } 
s<-svd(t(B[,ah]));v[ah,]<-s$u[,1:r]%*%t(s$v)
        bb<-bb+sum(s$d)
}

 cc<-(sqrt(aa)+sqrt(bb))/2
}

 for (k in 1:r) {
   for (i in 1:M) {
    ak<-(ck[i]+1):ck[i+1]
     for (j in 1:N) {
     ah<-(ch[j]+1):ch[j+1]
    s2[i,j,k]<-(t(u[ak,k])%*%x[ak,ah]%*%v[ah,k])^2
     }
   }
 }
list(u=u,v=v,s2=s2)
}

