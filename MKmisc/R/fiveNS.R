fiveNS <- function(x, na.rm = TRUE, type = 7){
  xna <- is.na(x)

  if (na.rm) 
    x <- x[!xna]
  else if (any(xna)) 
    return(rep.int(NA, 5))

  n <- length(x)
  if (n == 0) 
    return(rep.int(NA, 5))
  else{
    return(c('Minimum' = min(x, na.rm = na.rm), 
             quantile(x, prob = 0.25, type = type, na.rm = na.rm),
             'Median' = median(x, na.rm = na.rm),
             quantile(x, prob = 0.75, type = type, na.rm = na.rm),
             'Maximum' = max(x, na.rm = na.rm)))
  }
}
