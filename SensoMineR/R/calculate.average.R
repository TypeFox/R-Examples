"calculate.average" <- function(axeACP){
nbcoord <- ncol(axeACP)-2
nbjuge <- -1+length(levels(axeACP[,dim(axeACP)[[2]]]))
nbprod <- length(levels(axeACP[,dim(axeACP)[[2]]-1]))
  tab<-array(0,dim=c(nbjuge,nbprod,nbcoord))
  moy<-matrix(0,nbcoord,nbprod)
  for (j in 1:nbprod){
  for (k in 1:nbcoord){
    moy[k,j]<-axeACP[j,k]
    }}
  for (k in 1:nbcoord){
  for (i in 1:nbjuge){
  for (j in 1:nbprod){
    tab[i,j,k]<-axeACP[nbprod+((i-1)*nbprod+j),k]
    }}}
  calc.moy <- list()
  calc.moy$tab <- tab
  calc.moy$moy <- moy
  return(calc.moy)
}
