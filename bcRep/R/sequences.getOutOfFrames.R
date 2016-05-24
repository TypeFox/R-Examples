## Julia Bischof
## 10-09-2015

sequences.getOutOfFrames<-function(data=NULL){
  if(length(data)==0 || length(grep("JUNCTION_frame",colnames(data)))==0){
    stop("--> Data frame is missing or it contains no Junction frame information")
  }
  data<-data[grep("out-of-frame",data$JUNCTION_frame),]
  return(data)
}

