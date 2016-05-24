"fastflightcube" <-
function(X,pik,order=1,comment=TRUE) 
{ 
EPS = 1e-11
"algofastflightcube" <-
function(X,pik) 
{ 

"jump" <-
function(X,pik){ 
N = length(pik)
p = round(length(X)/length(pik))
X<-array(X,c(N,p))
X1=cbind(X,rep(0,times=N)) 
kern<-svd(X1)$u[,p+1] 
listek=abs(kern)>EPS
buff1<-(1-pik[listek])/kern[listek]
buff2<- -pik[listek]/kern[listek]
la1<-min( c(buff1[(buff1>0)] , buff2[(buff2>0)]) ) 
pik1<- pik+la1*kern 
buff1<- -(1-pik[listek])/kern[listek] 
buff2<- pik[listek]/kern[listek]
la2<-min(c(buff1[(buff1>0)] , buff2[(buff2>0)])) 
pik2<- pik-la2*kern 
q<-la2/(la1+la2)  
if (runif(1)<q) pikn<-pik1 else pikn<-pik2 
pikn 
}

N = length(pik)
p = round(length(X)/length(pik))
X<-array(X,c(N,p))
A<- X/pik 
B<-A[1:(p+1),]
psik <- pik[1:(p+1)]
ind<-seq(1,p+1,1) 
pp=p+2 
B<-array(B,c(p+1,p))
while(pp<=N) 
{ 
psik <- jump(B,psik)
liste<- (psik>(1-EPS) | psik<EPS)
i<- 0 
while(i <=(p) & pp<=N) 
{ i=i+1
if(liste[i]==TRUE) 
{pik[ind[i]]=psik[i]
psik[i]=pik[pp] 
B[i,]=A[pp,]
B=array(B,c(p+1,p))
ind[i]=pp 
pp=pp+1 } 
} 
} 
if(length(pik[(pik>EPS & pik<(1-EPS))])==(p+1)) psik <- jump(B,psik)
pik[ind]=psik
pik 
}
"reduc" <-
function(X)
{
EPS=0.0000000001
N=dim(X)[1]
Re=svd(X)
array(Re$u[,(Re$d>EPS)] , c(N,sum(as.integer(Re$d>EPS))))
}

N = length(pik);
p = round(length(X)/length(pik))
X<-array(X,c(N,p))
if (order==1) o<-sample(N,N)  else
   {
   if(order==2) o<-seq(1,N,1) 
    else o<-order(pik,decreasing=TRUE)
    }
liste<-o[(pik[o]>EPS & pik[o]<(1-EPS))]
if(comment==TRUE){
cat("\nBEGINNING OF THE FLIGHT PHASE\n")
cat("The matrix of balanced variable has",p," variables and ",N," units\n")
cat("The size of the inclusion probability vector is ",length(pik),"\n")
cat("The sum of the inclusion probability vector is ",sum(pik),"\n")
cat("The inclusion probability vector has ",length(liste)," non-integer elements\n")
} 
pikbon<-pik[liste]; 
Nbon=length(pikbon); 
Xbon<-array(X[liste,] ,c(Nbon,p)) 
pikstar<-pik 
flag=0
if(Nbon>p){if(comment==TRUE) cat("Step 1  ")
           pikstarbon<-algofastflightcube(Xbon,pikbon)
           pikstar[liste]=pikstarbon 
           flag=1
           }
liste<-o[(pikstar[o]>EPS & pikstar[o]<(1-EPS))]
pikbon<-pikstar[liste] 
Nbon=length(pikbon) 
Xbon<-array(X[liste,] ,c(Nbon,p))
pbon=dim(Xbon)[2]
if(Nbon>0){
          Xbon=reduc(Xbon)
          pbon=dim(Xbon)[2]
          }
k=2
while(Nbon>pbon & Nbon>0){
           if(comment==TRUE) cat("Step ",k,",  ")
           k=k+1
           pikstarbon<-algofastflightcube(Xbon/pik[liste]*pikbon,pikbon)
           pikstar[liste]=pikstarbon
           liste<-o[(pikstar[o]>EPS & pikstar[o]<(1-EPS))]
           pikbon<-pikstar[liste] 
           Nbon=length(pikbon) 
           Xbon<-array(X[liste,] ,c(Nbon,p))
           if(Nbon>0)
               {
               Xbon=reduc(Xbon)
               pbon=dim(Xbon)[2]
               }
           flag=1
           }
if(comment==TRUE) if(flag==0) cat("NO FLIGHT PHASE")
if(comment==TRUE) cat("\n")
pikstar 
}

