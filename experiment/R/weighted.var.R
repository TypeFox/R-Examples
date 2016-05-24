###
### weighted variance formula
###

weighted.var <- function(x, w) 
  return(sum(w * (x - weighted.mean(x,w))^2)/((length(x)-1)*mean(w)))
  
