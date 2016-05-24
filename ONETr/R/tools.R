tools <-
function(list){
  len <- 
  if(length(list$tools_technology$tools) > 0){
    tools <- ldply((lapply(list$tools_technology$tools, function(x){t(unlist(x))})))
    return(tools)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
