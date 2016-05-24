
validate.outcome <- function(null, resp.var){
  
  resp.level <- unique(null[, resp.var])
  
  if(length(resp.level) == 2){
    if(!setequal(resp.level, c(0, 1))){
      msg <- "response variable in formula should be 0 and/or 1"
      stop(msg)
    }
  }else{
    msg <- "response variable in formula should have two levels"
    stop(msg)
  }
  
}
