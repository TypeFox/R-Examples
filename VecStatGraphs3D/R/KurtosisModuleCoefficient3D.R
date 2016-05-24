KurtosisModuleCoefficient3D <- function (modules) 
{
  n = length(modules)  
  mean = sum(modules) / n
  s = ModuleStandardDeviation3D(modules)
  a = (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))
  b = sum(((modules - mean) / s)^4)
  cc = (3 * (n - 1)^2) / ((n - 2) * (n - 3))
  return((((n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))) * 
    (sum(((modules - mean) / s)^4))) 
    - ((3 * (n - 1)^2) / ((n - 2) * (n - 3))))
}
