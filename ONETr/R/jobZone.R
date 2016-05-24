jobZone <-
function(list){
  if(length(list$job_zone) > 0){
    job_zone <- ldply((lapply(list$job_zone, function(x){t(unlist(x))})))
    return(job_zone)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
