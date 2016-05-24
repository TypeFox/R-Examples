searchWordSD <- function(keeper, string, maxdist, summarize=F){
  bucket <- c()
  result <- list()
  buscar <- delList(string,maxdist)[[1]]
  for(item in buscar){
    if(searchTreeSD(keeper, item)[1]!=F && !is.null(searchTreeSD(keeper, item))){
      bucket<-append(bucket,searchTreeSD(keeper, item))
    }
  }
  bucket<-unique(bucket)
  bucket<-bucket[which(stringdist(string,bucket)<=maxdist)]
  if(length(bucket>0)){
    if(summarize==F){
      return(bucket)
    }
    else{
      if(length(bucket[which(stringdist(string,bucket)==1)])!=0){
        result$dist_1 <- bucket[which(stringdist(string,bucket)==1)]}
      if(length(bucket[which(stringdist(string,bucket)==2)])!=0){
        result$dist_2 <- bucket[which(stringdist(string,bucket)==2)]}
      if(length(bucket[which(stringdist(string,bucket)==3)])!=0){
        result$dist_3 <- bucket[which(stringdist(string,bucket)==3)]}
      return(result)
    }
  }
  else{NULL}
}