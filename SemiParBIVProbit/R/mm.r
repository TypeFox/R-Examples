mm <- function(ob){

  epsilon <- 0.0000001 
  max.p   <- 0.9999999

  res <- ifelse(ob > max.p, max.p, ob) 
  res <- ifelse(res < epsilon, epsilon, res) 

  res
      
}


