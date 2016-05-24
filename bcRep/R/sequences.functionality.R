## Julia Bischof
## 10-09-2015

sequences.functionality<-function(data=NULL, relativeValues=TRUE, ...){
  if(is.vector(data)==F || length(grep("productive",data))==0){
    stop("--> Functionality vector is missing or empty")
  }
  if(relativeValues==T){
    func.tab<-data.frame(length(grep("^productive",data))/length(data),
                         length(grep("^unproductive",data))/length(data),
                         (length(data)-length(grep("^productive",data))-length(grep("^unproductive",data)))/length(data))
  }else{
    func.tab<-data.frame(length(grep("^productive",data)),
                         length(grep("^unproductive",data)),
                         (length(data)-length(grep("^productive",data))-length(grep("^unproductive",data))))
  }
  colnames(func.tab)<-c("productive", "unproductive","unknown")
  return(data.frame(func.tab, check.names=F, row.names=NULL))
}

