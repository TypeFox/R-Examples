relatedOccupations <-
function(list){
  if(length(list$related_occupations) > 0){
    related_occupations <- ldply((lapply(list$related_occupations, function(x){t(unlist(x))})))
    return(related_occupations)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
