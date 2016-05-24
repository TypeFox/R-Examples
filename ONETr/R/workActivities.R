workActivities <-
function(list){
  if(length(list$work_activities) > 0){
    work_activities <- ldply((lapply(list$work_activities, function(x){t(unlist(x))})))
    return(work_activities)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
