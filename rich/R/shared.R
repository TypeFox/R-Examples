shared <-
function(lmatrix)
 {
data<-lmatrix
if(length(unique(as.vector(sapply(data,dim)[2,])))!=1) stop("some data matrices have different number of columns (species)")
n<-length(data)
nb.res<-n*(n-1)/2
matres<-matrix(data=NA,nrow=n,ncol=n)
  for (i in 1:(n)) {
      for (j in 1:(n)) {
      if(i==j) {
	  d<-data[[i]]
	  matres[i,j]<-length(colSums(d)[colSums(d)>0])}
      if(i>j){
	  a<-colSums(data[[i]])
	  b<-colSums(data[[j]])
	  am<-a;bm<-b
	  indice.des.doubles.0<-which(colSums(rbind(a,b))==0)
	  if(length(indice.des.doubles.0)>0){
	  am<-a[-indice.des.doubles.0];bm<-b[-indice.des.doubles.0]}
	  abin<-am;abin[abin>1]<-1;bbin<-bm;bbin[bbin>1]<-1
	  abin;bbin
	  shared.sp<-length(which(abin==bbin))
	  richesse.paire<-dim(data[[i]])[2]-length(indice.des.doubles.0)
	  matres[i,j]<-richesse.paire
	  matres[j,i]<-shared.sp}
      }
  }
rownames(matres)<-names(data)
colnames(matres)<-names(data)    
return(matres)    
}

