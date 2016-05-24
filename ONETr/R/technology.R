technology <-
function(list){
  if(length(list$tools_technology$technology) > 0){
    technology <- ldply((lapply(list$tools_technology$technology, function(x){t(unlist(x))})))
    return(technology)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
