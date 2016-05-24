abilities <-
function(list){
  if(length(list$abilities) > 0){
    abilities <- ldply((lapply(list$abilities, function(x){t(unlist(x))})))
    return(abilities)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
