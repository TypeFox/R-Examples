rcsn <- function(k, mu=rep(0, n), sigma, gamma, nu=rep(0, q), delta){
  if (!is.matrix(sigma) && length(mu) != 1)
    stop("incorrect 'mu'")
  if (is.matrix(gamma) && nrow(gamma) != length(nu))
    stop("incorrect 'gamma'")
  if (is.matrix(gamma) && ncol(gamma) != length(mu))
    stop("incorrect 'gamma'")
  if (is.matrix(sigma) && nrow(sigma) != length(mu))
    stop("incorrect 'sigma'")
  if (is.matrix(sigma) && ncol(sigma) != length(mu))
    stop("incorrect 'sigma'")
  if (is.matrix(delta) && nrow(delta) != length(nu))
    stop("incorrect 'delta'")
  if (is.matrix(delta) && ncol(delta) != length(nu))
    stop("incorrect 'delta'")
  n <- if (is.matrix(sigma)) 
    ncol(sigma)
  else 1
  res <- matrix(0,k,n)
  for (i in 1:k){
   repeat{
    v <- rmvnorm(1,mean = -nu, sigma= delta+gamma%*%sigma%*%t(gamma))
    if (all(v>=0)) break
  }
  temp <- t(gamma%*%sigma)%*%solve(delta+gamma%*%sigma%*%t(gamma))
  E <- mu+temp%*%t(v+nu)
  Var <- sigma-temp%*%(gamma%*%sigma)
  t <- rmvnorm(1, mean = E, sigma = as.matrix(Var))
  res[i,] <- t
  }
  return (res)
}
