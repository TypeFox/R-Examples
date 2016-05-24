SubsetterPIC <-
function(genos,save=NULL){
  genos<-as.matrix(genos)
  genos[is.na(genos)]<-0
  n.genos<-ncol(genos)
  result.list<-matrix(0,ncol=2,nrow=n.genos)
  colnames(result.list)<-c("Individual","PIC")
  result.list[1,2]<-PicCalc(genos)
  m<-nrow(genos)
  length.save<-length(save)
  if(length.save>1){
    stop<-length(save)+1
  }else{
    stop<-2
  }
  join<-array(0,c(2,m,n.genos))
  join[2,,]<-as.matrix(genos)
  mat.a<-colSums(join==1)
  mat.b<-colSums(join==-1)
  names<-colnames(genos)
  colnames(mat.a)<-names
  colnames(mat.b)<-names
  for(i in 1:(n.genos-stop)){
    names<-colnames(mat.a)
    a<-rowSums(mat.a)
    b<-rowSums(mat.b)
    temp.a<-a-mat.a
    temp.b<-b-mat.b
    temp.log<-colMeans(1-(temp.a^2+temp.b^2)/(temp.a+temp.b)^2,na.rm=TRUE)  #most intensive step- ~60% of time
  if(length(save)!=0){
  temp.log[save]<-0
  }
    remove<-which.max(temp.log)
    result.list[i,1]<-names[remove]
    result.list[i+1,2]<-max(temp.log)
    n.genos<-n.genos-1
    #removing is also about 30% of time.
    mat.a<-mat.a[,-which(colnames(mat.a) %in% names[remove])]
    mat.b<-mat.b[,-which(colnames(mat.b) %in% names[remove])]
  }
  length.genos<-nrow(result.list)
  result.list<-cbind(c(nrow(result.list):1),result.list)
  if(length.save>1){
    result.list[(length.genos-stop+1),2]<-setdiff(colnames(mat.a),save)
    result.list[(length.genos-stop+2):length.genos,2]<-save
    result.list[(length.genos-stop+2),3]<-PicCalc(genos[,save])
    result.list[(length.genos-stop+2):length.genos,1]<-1
  }else{
    result.list[(length.genos-1):length.genos,1]<-1
    result.list[(length.genos-1):length.genos,2]<-c(colnames(mat.a))
  }
  return(result.list)
}
