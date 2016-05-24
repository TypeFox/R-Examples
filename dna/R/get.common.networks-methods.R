setMethod("get.common.networks","pairOfNetworks",
function(object){
 N1=slot(object,"network1")
 N2=slot(object,"network2")
 if ((!is.null(colnames(N1)))&(!is.null(colnames(N2)))){
  common.genes=intersect(colnames(N1),colnames(N2))
  iN1=match(common.genes,colnames(N1))
  iN2=match(common.genes,colnames(N2))
  network1=N1[,iN1]
  network2=N2[,iN2]
 }
 else{
  network1=N1
  network2=N2
  p=ncol(N1)
  colnames(network1)=paste("Gene",1:p)
  colnames(network2)=paste("Gene",1:p)
 }
 return(list(network1=network1,network2=network2))
})

