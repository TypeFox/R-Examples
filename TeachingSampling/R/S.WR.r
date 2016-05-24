S.WR<-function(N,m){
nk<-rep(0,N)
  for(k in 1:N){
  suma<-sum(nk)
  nk[k]<-rbinom(1,(m-suma),(1/(N-k+1)))
               }
x<-which(nk>0)
w<-nk[x]
sam<-rep(x[1],w[1])

if(length(x)==1){
return(sam)}

if(length(x)>1){
for(i in 2:length(x)){
sam<-c(sam,rep(x[i],w[i]))
                     }
}
sam
}