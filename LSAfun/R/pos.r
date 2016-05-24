#######################################
##### set negative values to zero #####


pos <- function(x){
  
  x[x < 0] <- 0
  x
  
}