resampling_gauss <-
function(h, n, X)
{
  ind = sample(n, size = n, replace = TRUE)
  Z = rnorm(n, mean = 0, sd = 1)
  XB = h*Z + X[ind]
  
  XB = sort(XB)
  res <- XB
  return (res)
}
