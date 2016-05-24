relative_sample <-
function(X0, X1, h)
{

  n1 = length(X1)
  x <- numeric(length = n1)

  if (h == 0)
  {
      x = cdf_empirical_dist(X0,X1)
  }

  else
  {
      x = cdf_gaussianKernel(X0,X1,h)
  }

  res <- x
  return(res)
}
