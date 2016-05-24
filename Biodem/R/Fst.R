"Fst" <- function(rval,N){
  k<-N/sum(N)
  Fst.val<-k%*%diag(rval)
}
