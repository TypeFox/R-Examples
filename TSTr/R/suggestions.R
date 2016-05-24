suggestions<- function(tree,string, useUpper = F){
  vara <- variations(string, useUpper)
  bucket <- c()
  result <- list()
  for(item in vara){
    if(searchWord(tree,item) == T){
      bucket <- unique(append(bucket, item))
    }
  }
  return(bucket)
}
