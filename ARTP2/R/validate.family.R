
validate.family <- function(family){
  
  if(!(family %in% c('gaussian', 'binomial'))){
    msg <- 'family should be either \'gaussian\' or \'binomial\''
    stop(msg)
  }
  
}
