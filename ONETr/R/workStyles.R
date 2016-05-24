workStyles <-
function(list){
  if(length(list$work_styles) > 0){
    work_styles <- ldply((lapply(list$work_styles, function(x){t(unlist(x))})))
    return(work_styles)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
