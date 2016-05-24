
# this function is not efficient if data is huge
# will optimize it later
data.subset <- function(data, subset, options){
  
  if(!is.null(subset)){
    subset <- intersect(subset, 1:nrow(data))
    if(length(subset) > 0){
      data <- data[subset, ]
    }else{
      msg <- "No data left. Check subset"
      if(options$print) message(msg)
    }
  }
  
  data
  
}
