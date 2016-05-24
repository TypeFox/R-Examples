SubsetterMTK <-
function(genos,save=NULL,power=10,mat=NULL){
  genos<-as.matrix(genos)
  mode(genos)<-"numeric"
  genos[genos==0]<-NA
  n.genos<-ncol(genos)
  result.list<-matrix(0,ncol=4,nrow=n.genos)
  result.list[,1]<-n.genos:1
  colnames(result.list)<-c("Rank","Individual","Score","Mean Kinship")
  if(is.null(mat)){
    mat<-A.mat(t(genos))
    mat<-(2*(mat-min(mat))/(max(mat)-min(mat)))^power
    for(i in 1:n.genos){
      mat[i,i]<-0
    }
  }
  length.save<-length(save)
  if(length.save>1){
    stop<-length(save)
  }else{
    stop<-2
  }
  for(i in 1:(n.genos-stop)){
    result.list[i,3]<-sum(mat)/(n.genos^2)
    temp<-colSums(mat)
    if(length(save)!=0){
      ignore<-which(colnames(mat)%in%save)
      temp[ignore]<-0
    }
    remove<-which.max(temp)
    result.list[i,2]<-names(remove)
    mat<-mat[-remove,-remove]
  }
  if(length.save>1){
    result.list[(n.genos-stop+1):n.genos,2]<-save
    result.list[(n.genos-stop+1),3]<-sum(mat)/(n.genos^2)
    result.list[(n.genos-stop+1):n.genos,1]<-1
  }else{
    result.list[(n.genos-1):n.genos,1]<-1
    result.list[(n.genos-1):n.genos,2]<-c(colnames(mat))
    result.list[(n.genos-1),3]<-sum(mat)/(n.genos^2)
  }
   result.list[,4]<-as.numeric(result.list[,3])^(1/power)
  return(result.list)
}
