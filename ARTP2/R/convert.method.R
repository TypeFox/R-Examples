
convert.method <- function(method){
  
  if(length(method) > 1){
    msg <- 'Invalid options$method'
    stop(msg)
  }
  
  if(is.numeric(method)){
    return(method)
  }
  
  valid.method <- c('ADAJOINT', 'ADAJOINT2', 'ARTP')
  method <- toupper(method)
  if(!(method %in% valid.method)){
    msg <- 'Invalid options$method'
    stop(msg)
  }
  
  method <- which(valid.method == method)
  method
  
}
