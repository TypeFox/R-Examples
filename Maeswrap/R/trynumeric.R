trynumeric <- function(x){
  
  if(all(is.na(suppressWarnings(as.numeric(x)))))
    x
  else
    as.numeric(x)
  
}