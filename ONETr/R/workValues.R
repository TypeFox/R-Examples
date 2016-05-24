workValues <-
function(list){
  if(length(list$work_values) > 0){
    work_values <- ldply((lapply(list$work_values, function(x){t(unlist(x))})))
    return(work_values)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
