Combinations <-
function(n, k){
  # Compute all n choose k combinations of size k from 1:n
  # Return matrix with k rows and choose(n,k) columns.
  # Avoids recursion.  Code provided by Tim Hesterberg
  if(!is.numeric(n) || length(n) != 1 || n%%1) stop("'n' must be an integer")
  if(!is.numeric(k) || length(k) != 1 || k%%1) stop("'k' must be an integer")
  if(k > n || k <= 0) return(numeric(0))
  rowMatrix <- function(n) structure(1:n, dim=c(1,n))
  colMatrix <- function(n) structure(1:n, dim=c(n,1))
  if(k == n) return(colMatrix(n))
  if(k == 1) return(rowMatrix(n))
  L <- vector("list", k)
  # L[[j]] will contain combinations(N, j) for N = 2:n
  L[[1]] <- rowMatrix(2)
  L[[2]] <- colMatrix(2)
  Diff <- n-k
  for(N in seq(3, n, by=1)){
    # loop over j in reverse order, to avoid overwriting
    for(j in seq(min(k, N-1), max(2, N-Diff), by= -1))
      L[[j]] <- cbind(L[[j]], rbind(L[[j-1]], N, deparse.level=1))
    if(N <= Diff+1) L[[1]] <- rowMatrix(N)
    else L[[N-(Diff+1)]] <- numeric(0)
    if(N <= k) L[[N]] <- colMatrix(N)
  }
  L[[k]]
}

