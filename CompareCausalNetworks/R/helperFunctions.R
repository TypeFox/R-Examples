# function to determine number of settings to draw
drawE <- function(x){
  z <- floor(x)
  if(runif(1) < x-z) z <- z+1
  z <- max(1,z)
  return(z)
}

