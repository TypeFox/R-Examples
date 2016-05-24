# convert list of org.apache.mahout.cf.taste.recommender.RecommendedItem to an R list
.org.apache.mahout.cf.taste.recommender.RecommendedItem_list_to_Rlist <- function(x){  
  listsize <- x$size()  
  result <- list()
  result$items <- numeric(length = listsize)
  result$values <- numeric(length = listsize)
  i <- 0L
  while(i < listsize){
    entry <- x$get(i)    
    result$items[i+1] <- entry$getItemID()
    result$values[i+1] <- entry$getValue()
    i <- i+1L
  }
  result
}
