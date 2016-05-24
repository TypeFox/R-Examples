size.prop.test <- function(p=NULL,alpha,delta){
  u <- qnorm(1-alpha/2)^2
  if (is.null(p)) 
    {n <- u^2/(4*delta^2)} 
  else {n <- p*(1-p)*u^2/(delta^2)}
  return(ceiling(n))
}
