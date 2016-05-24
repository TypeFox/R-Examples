tap <- function(dataSet,one=NULL){
  
  # if test all pairs against one then extract the one
  if(!is.null(one)){
    z <- dataSet[,one]
    dataSet <- dataSet[,!(names(dataSet) %in% c(one))]
  }
  
  names <- names(dataSet)
  lr<-length(names)
  outMat = array(-1, c(lr,lr))
  print(paste("Processing",toString(lr),"variables.",toString(lr*(lr-1)/2),"pairs"))
  for (i in 1:(lr-1)) {
    for (j in (i+1):lr) {
      if(is.null(one)){
        d <- data.frame(V1=dataSet[,i],V2=dataSet[,j])
      }else{
        d <- data.frame(V1=dataSet[,i],V2=dataSet[,j],V3=z)
      }     
      dok <- d[complete.cases(d),]
      if(length(dok[,1])>10){
        aScore <- round(suppressWarnings(ma(dok)$A),2)
        outMat[i,j]<-aScore
        outMat[j,i]<-aScore  
      } else {
        outMat[i,j]<-0
        outMat[j,i]<-0
      }
    }
    print(paste(toString(i)," variable(s) complete",sep = ""))
  }  
  dFrame<-data.frame(outMat,row.names = names)
  names(dFrame)<-names
  return(dFrame)
}
