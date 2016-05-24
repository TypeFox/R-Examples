tdist <- function(df)  {
  ### df could be a vector or a scalar   
  
  ans <- list()   
  class(ans) <- "newfam"
  ans$family <- "tdist"
  ans$df <- df
  return(ans)
} 
