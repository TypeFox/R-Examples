CAP<-function(x, transform=NULL) {
  cap<-function(x) {
    y = as.data.frame(x)
    if(ncol(y)>1) {
      for(i in (ncol(y)-1):1) {
        y[,i] = y[,i]+y[,i+1]
      }
    }
    return(y)
  }
  if(!inherits(x,"stratifiedvegdata")) stop("Input should be of class 'stratifiedvegdata'")
  Y = lapply(x, FUN=cap)
  if(!is.null(transform)) Y = lapply(Y, FUN=transform)
  class(Y)<-c("list","CAP")
  return(Y)
}