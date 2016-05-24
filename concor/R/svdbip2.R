"svdbip2" <-
function(x,K,H,r) {
# INITIALISATIONS 
p<-dim(x)[1]
q<-dim(x)[2]

if (sum(H) != q | sum(K) != p) stop("K or H IS NOT SUITABLE")
if (r > min(c(K,H))) stop("r IS NOT SUITABLE")
M<-length(K)
N<-length(H)
u<-matrix(0,p,r);u1<-u
v<-matrix(0,q,r);v1<-v
ck<-cumsum(c(0,K))
ch<-cumsum(c(0,H))
s2<-array(0,c(M,N,r))

# INITIALISATIONS 

# ALGORITHM
for (k in 1:r) {
 a<-2
 #comp<-0

# PROPOSED INITIALISATION OF THE ALGORITHM, for u and v
 for (i in 1:M) {
   for (j in 1:N) {
    ak<-(ck[i]+1):ck[i+1]
    ah<-(ch[j]+1):ch[j+1]
    s<-svd(matrix(x[ak,ah],nrow=length(ak)))
    u[ak,k]<-s$u[,1];v[ah,k]<-s$v[,1]
   }
 }
a<-2;b<-0
 while (abs(a-b) > 1e-8) {
 #a^2 converge to the optimized criterion
 b<-a
  #comp<-comp+1
   a<-0
    v2<-matrix(0,q,1)
   for (j in 1:N) {
     ah<-(ch[j]+1):ch[j+1]
      for (i in 1:M) {
       ak<-(ck[i]+1):ck[i+1]
       v2[ah]<-v2[ah]+ t(matrix(x[ak,ah],nrow=length(ak)))%*%(u[ak,k]%*%t(u[ak,k])%*%x[ak,ah]%*%v[ah,k])
      }
     a2<-sqrt(t(v2[ah])%*%v2[ah]);
     if (a2 > 1e-8) v[ah,k]<-v2[ah]/a2 else v[ah,k]<-v2[ah]
   }

    u2<-matrix(0,p,1)
   for (i in 1:M) {
       ak<-(ck[i]+1):ck[i+1]
       for (j in 1:N) {
         ah<-(ch[j]+1):ch[j+1]
        u2[ak]=u2[ak]+ matrix(x[ak,ah],nrow=length(ak))%*%v[ah,k]%*%t(v[ah,k])%*%t(matrix(x[ak,ah],nrow=length(ak)))%*%u[ak,k]
       }
     a2<-sqrt(t(u2[ak])%*%u2[ak]);a<-a+a2
     if (a2 > 1e-8) u[ak,k]<-u2[ak]/a2  else u[ak,k]<-u2[ak]
   }

   a<-sqrt(a)
 }
 
   for (i in 1:M) {
    ak<-(ck[i]+1):ck[i+1]
     for (j in 1:N) {
     ah<-(ch[j]+1):ch[j+1]
     c<-t(u[ak,k])%*%x[ak,ah]%*%v[ah,k]    
    
x[ak,ah]<-x[ak,ah]-u[ak,k]%*%t(u[ak,k])%*%x[ak,ah]-x[ak,ah]%*%(v[ah,k]%*%t(v[ah,k]))+u[ak,k]%*%c%*%t(v[ah,k])   
     s2[i,j,k]<-c^2
     }
   }
 #comp
}
list(u=u,v=v,s2=s2)
}

