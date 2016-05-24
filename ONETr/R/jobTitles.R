jobTitles <- function(list){
  if(length(list$occupation$sample_of_reported_job_titles) > 0){
    jobTitles <- ldply((lapply(list$occupation$sample_of_reported_job_titles, function(x){t(unlist(x))})))
    return(jobTitles)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
