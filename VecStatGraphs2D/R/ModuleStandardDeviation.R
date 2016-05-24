ModuleStandardDeviation <- function (modules) 
{
  n = length(modules)
  mean = sum(modules) / n
  return(sqrt((sum((modules - mean)^2))/(n - 1)))
}
