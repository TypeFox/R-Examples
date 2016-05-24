## Julia Bischof
## 10-09-2015

sequences.junctionFrame<-function(data=NULL, relativeValues=TRUE, ...){
  if(is.vector(data)==F || length(grep("-frame",data))==0){
    stop("--> Junction frame vector is missing or empty")
  }
  if(relativeValues==T){
    junction.tab<-data.frame(length(grep("^in-frame",data))/length(data),
                             length(grep("^out-of-frame",data))/length(data),
                             (length(data)-length(grep("^in-frame",data))-length(grep("^out-of-frame",data)))/length(data))
  }else{
    junction.tab<-data.frame(length(grep("^in-frame",data)),
                             length(grep("^out-of-frame",data)),
                             (length(data)-length(grep("^in-frame",data))-length(grep("^out-of-frame",data))))
  }
  colnames(junction.tab)<-c("in-frame", "out-of-frame","unknown")
  return(data.frame(junction.tab, check.names=F, row.names=NULL))
}


