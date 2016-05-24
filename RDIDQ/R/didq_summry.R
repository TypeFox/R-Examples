didq_summry <-
function(test) {
  
  len=length(test)
  sm1<-matrix(1,nrow=len,ncol=6)
  nm=names(test)
  
  colnames(sm1)<-c("Variable Name","Type","No. of Distinct Values","missing Percentage","total Non Missing values","threshold")
  
  for (i in 1:ncol(test)){
    sm1[i,3]<-as.numeric(length(unique(test[,i])))
  }
  
  for( i in 1:ncol(test)){
    sm1[i,5]=as.numeric(sum(!is.na(test[,i])))
    sm1[i,6]=as.numeric(sm1[i,3]*100)/as.numeric(sm1[i,5])
  }
  
  #missing percentage
  nas <-colSums(is.na(test))
  p=(nas/nrow(test))*100
  sm1[,4]=unlist(p)
  mysum <- data.frame(sm1)
  mysum <- mysum
  #Assigning row names
  for (i in 1:ncol(test)){
    sm1[i,1]=nm[i]
  }
  
  #define the class
  cl <- lapply(test, class)
  sm1[,2]=unlist(cl)
  sum1<-as.data.frame(sm1)
  return(sum1[,-6])
  
}
