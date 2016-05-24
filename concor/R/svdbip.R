"svdbip" <-
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
A<-matrix(0,p,N)
B<-matrix(0,q,M)
s2<-array(0,c(M,N,r))

# ALGORITHM
for (k in 1:r) {

 #PROPOSED INITIALISATION OF THE ALGORITHM with u
    for (i in 1:M) {
      ak<-(ck[i]+1):ck[i+1]
     s<-svd(matrix(x[ak,],nrow=length(ak)))
      u[ak,k]<-s$u[,1]
    }


 cc<-s$d[1];cc1<-0;
 #comp<-0;

 while (abs(cc-cc1) > 1e-8) {
  #aa^2 and bb^2 are converging to the optimized criterion
  aa<-0;bb<-0;
  cc1<-cc;

  #comp<-comp+1;

   for (j in 1:N) {
    ah<-(ch[j]+1):ch[j+1]
     for (i in 1:M) {
      ak<-(ck[i]+1):ck[i+1]
      B[ah,i]<-t(matrix(x[ak,ah],nrow=length(ak)))%*%u[ak,k]
     }
    s<-svd(matrix(B[ah,],nrow=length(ah)))
    if (s$d[1] > 1e-8) { v[ah,k]<-s$u[,1]; aa<-aa+s$d[1] }
   }


   for (i in 1:M) {
     ak<-(ck[i]+1):ck[i+1]
       for (j in 1:N) {
        ah<-(ch[j]+1):ch[j+1]
      A[ak,j]<- matrix(x[ak,ah],nrow=length(ak))%*%v[ah,k]
       }
    s<-svd(matrix(A[ak,],nrow=length(ak)))
       if (s$d[1] > 1e-8) {u[ak,k]<-s$u[,1];bb<-bb+s$d[1]}
   }

  cc<-(aa+bb)/2
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

