## Julia Bischof
## 10-09-2015

sequences.getAnyFunctionality<-function(data=NULL){
  if(length(data)==0 || length(grep("Functionality",colnames(data)))==0){
    stop("--> Data frame is missing or it contains no Functionality information")
  }
  data<-data[grep("productive",data$Functionality),]
  return(data)
}

