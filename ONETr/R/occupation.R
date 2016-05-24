occupation <-
function(list){
  if(length(list$occupation) > 0){
    occupation <- as.data.frame(t(ldply((lapply(list$occupation[c(1,2,4)], function(x){t(unlist(x))})))))
    names(occupation) <- names(list$occupation[c(1,2,4)])
    occupation <- cbind(occupation, ldply((lapply(list$occupation[3], function(x){t(unlist(x))})))[,2:3])[2,]
    return(occupation)
  }
  else{
    message("Warning: This type of data is missing or incomplete for this occupation.")
  }
}
