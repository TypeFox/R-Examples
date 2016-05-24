logList <-  function(x){
    if(length(x)!=4){
      print("Please provide correct RPPanalyzer list format including 4 elements, i.e. matrix with protein expression data, matrix with background data, data frame with array features and data frame with sample description.")
    }
    x.log<-list()    
    x.log[[1]]<-log2(x[[1]])
    x.log[[2]]<-log2(x[[2]])
    x.log[[3]]<-x[[3]]
    x.log[[4]]<-x[[4]]
    return (x.log)
}