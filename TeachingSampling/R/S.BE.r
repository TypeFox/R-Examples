S.BE<-function(N,prob){
sam<-matrix(0,N,1)
U<-runif(N)
for(k in 1:N){
if(U[k]<=prob)
sam[k]<-k
              }
return(sam)
}