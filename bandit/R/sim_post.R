sim_post <- function(x,
 n,
 alpha = 1,
 beta = 1,
 ndraws = 5000
) {
  k <- length(x)
  ans <- matrix(nrow=ndraws, ncol=k)
  no = n-x
  for (i in (1:k))
    ans[,i] = rbeta(ndraws, x[i] + alpha, no[i] + beta)
  return(ans)
}
