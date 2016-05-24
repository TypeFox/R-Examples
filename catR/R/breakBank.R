breakBank<-function(itemBank){
nc<-ncol(itemBank)
itemPar<-matrix(as.numeric(as.matrix(itemBank[,1:(nc-1)])),nrow(itemBank),(nc-1))
cbGroup<-as.factor(itemBank[,nc])
res<-list(itemPar=itemPar,cbGroup=cbGroup)
return(res)}
