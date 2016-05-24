N2error <-
function(dat){
  nMarkers=dim(dat)[2]-1
  error=1
  names(error)=c("Max allele Count error")
  for (i in 1:nMarkers){
    p=dat[!is.na(dat[,i+1]),i+1]
    error=error*sum(N2Exact(p)[1:2])
  }
error
}

