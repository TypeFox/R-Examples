S.PPS<-function(m,x){
N<-length(x)
pk<-x/sum(x)
cumpk<-cumsum(pk)
U<-runif(m)
ints<-cbind(c(0,cumpk[-N]),cumpk)
sam<-rep(0,m)
for(i in 1:m){
    sam[i]<-which(U[i]>ints[,1] & U[i]<ints[,2])
   }
return(cbind(sam,pk[sam]))
              }