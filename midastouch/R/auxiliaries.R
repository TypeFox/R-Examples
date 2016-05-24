bootfunc.plain <- function(n){
  random <- sample(n,replace = TRUE)
  weights <- as.numeric(table(factor(random,levels = c(1:n))))
  return(weights)
}
minmax <- function(x,domin=TRUE,domax=TRUE){
  maxx <- sqrt(.Machine$double.xmax)
  minx <- sqrt(.Machine$double.eps)
  if(domin){      x <- pmin(x,maxx)   }
  if(domax){      x <- pmax(x,minx)   }
  return(x)
}
