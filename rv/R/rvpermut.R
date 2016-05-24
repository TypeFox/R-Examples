

# ========================================================================
# 
# ========================================================================

rvpermut <- function (data, prob=NULL) {
  ## permutation distribution
  n.sims <- getnsims()
  s <- t(sapply(rep(list(data), n.sims), sample, prob=prob))
  r <- rvsims(s)
  dim(r) <- dim(data)
  names(r) <- names(data)
  return(r)
}

