Amat <- function(type = "ridge", corrs, c = 1e-08){
  
  if(type == "lasso") A <- 10 # check
  if(type == "ridge") A <- diag(c(1, 1, 1))
   
A 
  
}