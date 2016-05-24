skills <-
function(list){
  if(length(list$skills) > 0){
    skills <- ldply((lapply(list$skills, function(x){t(unlist(x))})))
    return(skills)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
