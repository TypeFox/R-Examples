SupportRS <- function(N, ID=FALSE){
sam <- matrix(NA, ncol=N, nrow=1)
for(k in 1:N){
sam<-rbind(sam, cbind(Support(N,k),matrix(NA,ncol=N-k, nrow=choose(N,k))))
}
if(ID==FALSE){
return(sam)
}
else{
sam<-matrix(ID[SupportRS(N)],nrow=2^N)
return(sam)
}
}