SupportWR <- function(N, m, ID=FALSE){
S=0
a=rep(1,m)
P1<-a
S=S+1
k=m
while(k>0){
while(a[k]<N){
a[k]=a[k]+1
P1<-rbind(P1,a)
S=S+1
}
if(k>1)
k=k-1
if(a[k]<N){
a[k]=a[k]+1
k1=k+1
a[k1:m]=a[k]
P1<-rbind(P1,a)
S=S+1
k=m
}
else
if(k==1)
k=0
}

nr <- choose(N+m-1,m)
P1 <- matrix(P1, nrow=nr)
sam <- matrix(ID[P1], nrow=nr)

if (ID==FALSE) {return(P1)}
else {return(sam)}
}
