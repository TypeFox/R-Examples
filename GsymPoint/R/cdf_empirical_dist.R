cdf_empirical_dist <-
function(X, gridx)
{
  m = length(gridx)
  F <- numeric(length = m)
  for (i in 1:m)
  {
      F[i] = empirical_dist(X, gridx[i])
  }
  res <- F
  return(res)
}
