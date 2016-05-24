## Julia Bischof
## 10-09-2015

sequences.getUnproductives<-function(data=NULL){
  if(length(data)==0 || length(grep("Functionality",colnames(data)))==0){
    stop("--> Data frame is missing or it contains no Functionality information")
  }
  data<-data[grep("^unproductive",data$Functionality),]
  return(data)
}

