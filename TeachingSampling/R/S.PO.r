S.PO<-function(N,Pik){
sam<-matrix(0,N,1)
U<-runif(N)
for(k in 1:N){
if(U[k]<=Pik[k])
sam[k]<-k
              }
return(sam)
}