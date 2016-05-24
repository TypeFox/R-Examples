`LoglikelihoodAR` <-
function(phi, z, MeanValue=0){
if(length(phi)==0) phi=0
phis<-c(1,-phi)
y<-z-MeanValue
n<-length(z)
-log(DetAR(phi))/2 - (n/2)*log(sum(crossprod(phis,ChampernowneD(y,length(phis)-1,MeanZero=TRUE))*phis)/n)
}

