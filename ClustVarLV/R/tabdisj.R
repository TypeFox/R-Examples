tabdisj<-function (dataquali){
  
  dataquali<-as.data.frame(dataquali)
  n<-nrow(dataquali)
  p<-ncol(dataquali)
  matquali<-c()
  for(j in 1:p){
    fact<-as.factor(dataquali[,j])
    tab <- matrix(0, n, nlevels(fact))
    rownames(tab)<-rownames(dataquali)
    for (niv in 1:nlevels(fact)){
      tab[which(fact==levels(fact)[niv]),niv]<-1
    }
    colnames(tab)<-paste(colnames(dataquali)[j],levels(fact),sep=".")
    matquali=cbind(matquali,tab) 
  }
  return(matquali)

 }
 