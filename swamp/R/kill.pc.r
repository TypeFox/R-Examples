"kill.pc"<-function(g,pc,imputeknn=F,center=T){
  if(class(g)!="matrix"){stop("g is not a matrix")}
  if(!is.numeric(pc)){stop("pc is not numeric")}
  if(any(pc>ncol(g))){stop("one element of pc is larger than ncol(g)")}
  if(any(pc>nrow(g))){stop("one element of pc is larger than nrow(g)")}
  
  if (imputeknn==T){
         require(impute)
         gimpute<-impute.knn(g)
         g<-gimpute$data
         }
  
  if(center==T){s<-svd(g-rowMeans(g))}
  if(center==F){s<-svd(g)}
  D <- diag(s$d) 
  for (i in 1:length(pc)){
        D[pc[i],pc[i]]<-0 }

  g3<-s$u %*% D %*% t(s$v)
  
  dimnames(g3)<-dimnames(g)
  return(g3)
 }