ELquantiles_gsym <-
function(y, q0, X0, X1, n0, n1, h0, h1, rho)
{
  q1 = rho*(1-q0)
  w0 = 1/n0
  h0_inv = 1/h0
  dif0 = (y-X0)*h0_inv
  p0 = sum(kernel_cdf(dif0, "gaussian")*w0)

  if (q0 > 0 & q0 < 1 & p0 > 0 & p0 < 1)
  {
    loglik0 = 2*n0*(p0*log(p0/q0)+(1-p0)*log((1-p0)/(1-q0)))
  }
  else if (p0 == 0 & (q0 > 0 & q0 < 1))
  {
    loglik0 = 2*n0*((1-p0)*log((1-p0)/(1-q0)))
  }
  else if (p0 == 1 & (q0 > 0 & q0 < 1))
  {
    loglik0 = 2*n0*(p0*log(p0/q0))
  }
  else
  {
    loglik0 = 2*n0*1000
  }

  w1 = 1/n1
  h1_inv = 1/h1
  dif1 = (y-X1)*h1_inv

  p1 = sum(kernel_cdf(dif1, "gaussian")*w1)

  if (q1 > 0 & q1 < 1 & p1 > 0 & p1 < 1)
  {
    loglik1 = 2*n1*(p1*log(p1/q1)+(1-p1)*log((1-p1)/(1-q1)))
  }
  else if (p1 == 0 & (q1 > 0 & q1 < 1))
  {
    loglik1 = 2*n1*((1-p1)*log((1-p1)/(1-q1)))
  }
  else if (p1 == 1 & (q1 > 0 & q1 < 1))
  {
    loglik1 = 2*n1*(p1*log(p1/q1))
  }
  else
  {
    loglik1 = 2*n1*1000
  }

  res <- loglik0 + loglik1
  return(res)
}
