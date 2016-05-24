knowledge <-
function(list){
  if(length(list$knowledge) > 0){
    knowledge <- ldply((lapply(list$knowledge, function(x){t(unlist(x))})))
    return(knowledge)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
