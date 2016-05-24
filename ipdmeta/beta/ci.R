ci <- function(C,estimates,variance,alpha=.05,f=function(x){exp(x)}){

z <- qnorm(1-alpha/2)

g <- function(x){

  point.est = f(t(x)%*%estimates)
  low = f(t(x)%*%estimates-z*sqrt(t(x)%*%variance%*%x))
  high = f(t(x)%*%estimates+z*sqrt(t(x)%*%variance%*%x)) 
  
  result = array(c(low,point.est,high))
  names(result) = c("low","point.est","high")
  return(result)
  }

if(!is.matrix(C)){
  g(C)
 }
else{
  apply(C,1,g)
  }
}
