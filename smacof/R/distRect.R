distRect <- function(x,y,reg)
{
  return(sqrt(reg+outer(rowSums(x^2),rowSums(y^2),"+")-2*x%*%t(y)))
}
