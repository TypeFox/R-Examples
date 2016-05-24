tasks <-
function(list){
  if(length(list$tasks) > 0){
    tasks <- ldply((lapply(list$tasks, function(x){t(unlist(x))})))
    return(tasks)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
