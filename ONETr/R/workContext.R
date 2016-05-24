workContext <-
function(list){
  if(length(list$work_context) > 0){
    work_context <- ldply((lapply(list$work_context, function(x){t(unlist(x))})))
    return(work_context)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
