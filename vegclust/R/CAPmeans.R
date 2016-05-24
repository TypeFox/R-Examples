CAPmeans<-function(CAP, factor=NULL) {
  averageCAP<-function(x) {
    avccf = x[[1]]
    if(length(x)>1) {
      for(i in 2:length(x)) {
        avccf = avccf + x[[i]]
      }
      avccf = avccf/length(x)
    }
    return(avccf)
  }
  if(!is.null(factor)) {
    Y = lapply(split(CAP, factor), FUN = averageCAP)
  } else {
    Y = list(meanCAP=averageCAP(CAP))
  }
  class(Y)<-c("list","CAP")
  return(Y)  
}