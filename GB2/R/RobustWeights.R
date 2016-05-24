robwts <- function(x, w=rep(1, length(x)), c=0.01, alpha=0.001) {
  fiskpar <- fisk(x,w)  
  a <- fiskpar[1]                    
  b <- fiskpar[2]  
  num <- abs(((1-alpha)/alpha)^(1/a)-(alpha/(1-alpha))^(1/a))
  corr <- pmax(c,pmin(1,num/abs(b/x-1),num/abs(x/b-1)))
  # a list containing the correction and the adjusted weights
  return(list(corr, corr*w)) 
}

# prob <- 0.001 ;  c <- 0.1
