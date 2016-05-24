nchar0 <- function(x, ...){
  is.null(x) || (length(x) == 0) || 
    (max(nchar(x)) == 0) 
  
}
