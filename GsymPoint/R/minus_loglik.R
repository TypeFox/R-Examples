minus_loglik <-
function(lambda, parameters)
{
  n0 = parameters[1]
  X0 = parameters[2:c(n0+1)]
  N = length(lambda)
  n1 = parameters[c(n0+2)]
  X1 = parameters[c(n0+3):c(n0+n1+2)]
  g <- numeric(length = N)
  for (k in 1:N)
  {
    Y0 = ((X0^(lambda[k]))-1)/lambda[k]
    var0 = var(Y0)*(n0-1)/n0
    Y1 = ((X1^(lambda[k]))-1)/lambda[k]
    var1 = var(Y1)*(n1-1)/n1
    g[k] = 0.5*n0*log(2*pi)+0.5*n0*log(var0)-(lambda[k]-1)*(sum(log(X0)))+0.5*n1*log(2*pi)+0.5*n1*log(var1)-(lambda[k]-1)*(sum(log(X1)))
  }
  res <- g
  return(g)
}
