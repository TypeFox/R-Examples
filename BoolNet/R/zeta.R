
# Distribution function of the Zeta distribution
dzeta <- function(k,gamma=2.5,approx_cutoff=100)
{
  zeta_approx <- sum(sapply(seq_len(approx_cutoff),function(x)1/(x^gamma)))
  sapply(k,function(k_i)
  {
    1/(k_i^gamma * zeta_approx)
  })
}

# Quantile function of the Zeta distribution
qzeta <- function(p,maxK,gamma=2.5,approx_cutoff=100)
{
  vals <- cumsum(dzeta(seq_len(maxK),gamma,approx_cutoff))
  sapply(p,function(p_i)
  {
    indices <- which(vals > p_i)
    if (length(indices) == 0)
      return(NA)
    else
      return(indices[1])
  })
}

# Draw random numbers from the Zeta distribution
rzeta <- function(n,maxK,gamma=2.5,approx_cutoff=100)
{
  maxP <- cumsum(dzeta(seq_len(maxK),gamma,approx_cutoff))
  p <- runif(min=0,max=maxP,n=n)
  return(qzeta(p,maxK,gamma,approx_cutoff))
}
