`FastLoglikelihoodAR` <-
function(phi, n, CD){
phis<-c(1,-phi)
LL <- -log(DetAR(phi))/2-(n/2)*log(sum(crossprod(phis,CD)*phis)/n)
if (!is.finite(LL)){
   LL<--1e35 }
LL
}

