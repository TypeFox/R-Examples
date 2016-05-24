education <-
function(list){
  if(length(list$education$level_required) > 0){
    education <- ldply((lapply(list$education$level_required, function(x){t(unlist(x))})))
    return(education)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
