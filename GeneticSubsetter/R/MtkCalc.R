MtkCalc <-
function(genos,subset,power=10){
  genos<-as.matrix(genos)
  mode(genos)<-"numeric"
  genos[genos==0]<-NA
  n.genos<-ncol(genos)
  mat<-A.mat(t(genos))
  mat<-(2*(mat-min(mat))/(max(mat)-min(mat)))^power
  for(i in 1:n.genos){
    mat[i,i]<-0
  }
  return(sum(mat[subset,subset])/(n.genos^2))
}
