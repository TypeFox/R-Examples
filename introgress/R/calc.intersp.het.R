calc.intersp.het <-
function(introgress.data=NULL){
  if (is.null(introgress.data)==TRUE)
    stop("error, input data were not provided")
  if (is.list(introgress.data)==TRUE){
    count.matrix<-introgress.data[[2]]
  }
  else count.matrix<-introgress.data
  n.ind<-dim(count.matrix)[2]
  int.het<-numeric(n.ind)
  for (j in 1:n.ind){
    int.het[j]<-sum(count.matrix[,j]==1,na.rm=TRUE)/sum(is.na(count.matrix[,j])==FALSE)
  }	
  return(int.het)
}

