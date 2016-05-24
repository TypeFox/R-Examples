MTK<-
function(genos,subset,mat=NULL,power=10){
  genos<-as.matrix(genos)
  mode(genos)<-"numeric"
  n.genos<-ncol(genos)
  #Make and transform kinship matrix
  if(is.null(mat)){
    mat<-A.mat(t(genos))
  }
  mat.adj<-(2*(mat-min(mat))/(max(mat)-min(mat)))^power
  for(i in 1:n.genos){
    mat.adj[i,i]<-0
  }
  #Return mean of kinship matrix
  return(sum(mat.adj[subset,subset])/(n.genos^2))
}
