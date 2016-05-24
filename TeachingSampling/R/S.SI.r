S.SI<-function(N,n,e=runif(N))
{
c<-matrix(0,N,1)
dec<-matrix(0,N,1)
sam<-matrix(0,N,1)
for(k in 1:N){
  c[k]<-(n-dec[k])/(N-k+1)
  if(e[k]<c[k]){
  dec[k:N]<-dec[k]+1
  sam[k]<-k}
  }
sam
}