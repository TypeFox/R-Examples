interests <-
function(list){
  if(length(list$interests) > 0){
    interests <- ldply((lapply(list$interests, function(x){t(unlist(x))})))
    return(interests)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}