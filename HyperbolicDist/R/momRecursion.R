### Calculate recursion for any generalized hyperbolic distribution
### Christine Yang Dong and David Scott from code by Diethelm Wuertz
momRecursion <- function(order = 12, printMatrix = FALSE) {
   ## Description:
   ##   Computes the moment coefficients recursively
   ## Setting Start Values:
  a <- matrix(rep(0, times = order*order), ncol = order)
  a[1, 1] <- 1
  if (order > 1) {
    a[2, 1] <- 1
  }
  ## Compute all Cofficients by Recursion:
  if (order > 1) {
    for (d in 2:order) {
      for (l in 2:d) {
        a[d,l] <- a[d - 1,l - 1] + a[d - 1, l]*(2*l + 1 - d)
      }
    }
  }
  rownames(a) <- paste("order=", 1:order, sep = "")
  colnames(a) <- paste("l=", 1:order, sep = "")
  ## Print the matrix:
  if (printMatrix) {
    cat("\n")
    print(a)
    cat("\n")
  }
  for (k in 1:order) {
    L <- trunc((k + 1)/2):k
    M <- 2*L - k
  }
  return(list(a = a[order, L], L = L, M = M,
              lmin = trunc((order + 1)/2)))
}
