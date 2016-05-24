BerryMielke <- function(V) {
  # IOV, as cited in Blair & Lacy 2000
  # argument: V = frequency vector
    N <- sum(V)    # sample size
    k <- length(V) # number of categories
  # t  = sum from i=1 to k-1 ( sum from j = i+1 to k (ni*nj*(i-j)) )
  # ni, nj = cell counts; i, j = values of categories
    t <- sum(sapply(1:(k-1), function(i) sum(sapply((i+1):k, function(j) V[i]*V[j]*(j-i)))))
  # Tmax not defined with odd numbers; use (N^2-1)*(k-1)/4 instead (although differences small)
    Tmax <- N^2*(k-1)/4
  # if N odd:
    if(N%%2 == 1) (Tmax <- (N^2-1)*(k-1)/4) # %%: integer division, gives 1 if odd, gives 0 if even
  iov <- t/Tmax
  return(iov)
}
