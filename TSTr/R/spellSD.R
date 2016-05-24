spellSD<-function(keeper, string, maxdist=1, summarize=F){
  V2 <- c()
  buscar <- c()
  result <- c()
  bucket <- list()
  if(maxdist==1){
    buscar <- deletion1(string)
  }
  else if(maxdist==2){
    buscar <- deletion2(string)
  }
  else if(maxdist==3){
    buscar <- deletion3(string)
  }
  else{stop("argument \"maxdist\" out of range (1:3)")}
  for(item in buscar){
    result <- na.omit(unique(append(result,keeper[item,V2])))
  }
  if(!is.na(result[1])){
    if(summarize==F){
      result <- result[which(stringdist(string,result)<=maxdist)]
      return(result)
    }
    else{
      if(length(result[which(stringdist(string,result)==1)])!=0){
        bucket$dist_1 <- result[which(stringdist(string,result)==1)]}
      if(length(result[which(stringdist(string,result)==2)])!=0){
        bucket$dist_2 <- result[which(stringdist(string,result)==2)]}
      if(length(result[which(stringdist(string,result)==3)])!=0){
        bucket$dist_3 <- result[which(stringdist(string,result)==3)]}
      return(bucket)
    }
  }
  else{NULL}
}