cdf_gaussianKernel <-
function(X,gridx,h)
{
  n0 = length(X)
  m = length(gridx)
  w = 1/n0
  F <- numeric(length = m)
  dif <- numeric(length = n0)
  hinv = 1/h
  for (i in 1:m)
  {
    dif = (gridx[i]-X)*hinv
    F[i] = sum(kernel_cdf(dif,"gaussian")*w)
  }
  res <- F
  return (res)
}
