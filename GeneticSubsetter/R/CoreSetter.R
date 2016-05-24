CoreSetter<-
function(genos,criterion=c("HET","MTK"),save=NULL,power=10,mat=NULL){
  genos<-as.matrix(genos)
  mode(genos)<-"numeric"
  n.genos<-ncol(genos)
  length.save<-length(save)
  if(length.save>1){
    stop<-length(save)
  }else{
    stop<-2
  }
  ###Subsetting using HET criterion
  if(criterion=="HET"){
    result.list<-matrix(0,ncol=3,nrow=n.genos)
    colnames(result.list)<-c("Rank","Individual","HET")
    result.list[,1]<-c(n.genos:1)
    result.list[1,3]<-HET(genos)
    m<-nrow(genos)
    ##Build two matricies for efficient processing
    mat.a<-matrix(0,nrow=m,ncol=n.genos)
    mat.b<-matrix(0,nrow=m,ncol=n.genos)
    mat.a[genos==1]<-1
    mat.a[genos==0]<-0.5
    mat.b[genos==-1]<-1
    mat.b[genos==0]<-0.5
    names<-colnames(genos)
    colnames(mat.a)<-names
    colnames(mat.b)<-names
    ##Systematically remove genotypes
    for(i in 1:(n.genos-stop)){
      names<-colnames(mat.a)
      a<-rowSums(mat.a)
      b<-rowSums(mat.b)
      temp.a<-a-mat.a
      temp.b<-b-mat.b
      temp.log<-colMeans(1-(temp.a^2+temp.b^2)/(temp.a+temp.b)^2,na.rm=TRUE)
      if(length(save)!=0){
        temp.log[save]<-0
      }
      remove<-which.max(temp.log)
      result.list[i,2]<-names[remove]
      result.list[i+1,3]<-max(temp.log)
      mat.a<-mat.a[,-which(colnames(mat.a) %in% names[remove])]
      mat.b<-mat.b[,-which(colnames(mat.b) %in% names[remove])]
    }
    ##Clean Up
    if(length.save>1){
      result.list[(n.genos-stop+1):n.genos,1]<-1
      result.list[(n.genos-stop+1):n.genos,2]<-save
    }else{
      result.list[(n.genos-1):n.genos,1]<-1
      result.list[(n.genos-1):n.genos,2]<-c(colnames(mat.a))
    }
  }
  ###Subsetting using MTK Criterion
  if(criterion=="MTK"){
    result.list<-matrix(0,ncol=3,nrow=n.genos)
    result.list[,1]<-n.genos:1
    colnames(result.list)<-c("Rank","Individual","Score")
    if(is.null(mat)){
      mat<-A.mat(t(genos)) #Building Kinship Matrix
    }
    mat.adj<-(2*(mat-min(mat))/(max(mat)-min(mat)))^power #Transforming Kinship Matrix to Emphasize Similar Genotypes
    for(i in 1:n.genos){
      mat.adj[i,i]<-0
    }
    ##Systematically Remove Genotypes
    for(i in 1:(n.genos-stop)){
      result.list[i,3]<-sum(mat.adj)/(n.genos^2)
      temp<-colSums(mat.adj)
      if(length(save)!=0){
        ignore<-which(colnames(mat.adj)%in%save)
        temp[ignore]<-0
      }
      remove<-which.max(temp)
      result.list[i,2]<-names(remove)
      mat.adj<-mat.adj[-remove,-remove]
    }
    ##Clean Up
    if(length.save>1){
      result.list[(n.genos-stop+1):n.genos,2]<-save
      result.list[(n.genos-stop+1),3]<-sum(mat.adj)/(n.genos^2)
      result.list[(n.genos-stop+1):n.genos,1]<-1
    }else{
      result.list[(n.genos-1):n.genos,1]<-1
      result.list[(n.genos-1):n.genos,2]<-c(colnames(mat.adj))
      result.list[(n.genos-1),3]<-sum(mat.adj)/(n.genos^2)
    }
  }
  return(result.list)
}
