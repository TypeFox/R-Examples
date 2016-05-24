mutualCoherence <- function(A, k) {
  # Part of R1Magic by mehmet.suzen@physics.org
  # R recommenden library utils needed for combn
  # Single inner product
  .oneMutual <- function(Ap, Aq) {
    return((norm(t(Ap) %*% Aq, type='2')/(norm(Ap, type='2') * norm(Aq, type='2'))))
  }

  # Inner produce from two colums of a matrix
  .oneMutualcolIndex <- function(A, p, q) {
    return(.oneMutual(A[,p], A[,q]))
  }

  # Generate wrapper function on oneMutualcolIndex for known A
  .oneMutualcolIndexf <- function(A) {
    return(function(p, q) { .oneMutualcolIndex(A, p, q) } )
  }

  # Given set of colums ids Q, find the sum of the inner products
  .setOneMutual <- function(A, p, setQi) {
    mutualF <- .oneMutualcolIndexf(A)
    setQi   <- setdiff(setQi, p)
    pp      <- rep(p, length(setQi))
    vecQ    <- c()
    if(length(setQi)>0) vecQ    <- mapply(mutualF, pp, setQi)
    return(sum(vecQ))
  }

  # Generate wrapper function on setOneMutual for known A and p
  .setOneMutualf <- function(A, p) {
    return(function(setQi) { .setOneMutual(A, p, setQi) })
  }

  # Value of the max inner product in Qs; for a given column
  .setOneMutualP <- function(A, p, k) {
    nColsA     <- ncol(A)
    setQ       <- t(combn(nColsA, k)) # matlab - nchoosek
    oneMutualf <- .setOneMutualf(A, p)
    sumsQ      <- apply(setQ, 1, oneMutualf)
    return(max(sumsQ))
  }

  # Generate wrapper function on setOneMutalP for known A and k
  .setOneManualPf <- function(A, k) {
    return(function(p) { .setOneMutualP(A, p, k) } )
   }

  # Compute mutual coherence for Q sub-set of size k
  .kMutual <- function(A, k) {
     setOneMutualPff <- .setOneManualPf(A, k)
     nColsA          <- ncol(A)
     return(max(mapply(setOneMutualPff, 1:nColsA)))
  }

  # Generate wrapper function on kMutual for known A
  .kMutualf <- function(A, k) {
    return(function(k) { .kMutual(A, k) })
  }

   if(k > ncol(A)) stop("Invalid k, larger then number of columns")
   kMutualff <- .kMutualf(A, k)
   cohVec    <- mapply(kMutualff, 1:k)
   return(cohVec);
}
