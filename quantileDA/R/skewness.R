skewness<-function (x) 
{
  if (is.matrix(x)) 
    apply(x, 2, skewness)
  else if (is.vector(x)) {
    
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  }
  else if (is.data.frame(x)) 
    sapply(x, skewness)
  else skewness(as.vector(x))
}