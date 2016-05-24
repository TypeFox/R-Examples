# Computes x = A^+ * b using an SVD
# (slow but stable)

svdsolve <- function(A,b,rtol) {
  s = svd(A)
  di = s$d
  ii = di>rtol
  di[ii] = 1/di[ii]
  di[!ii] = 0
  return(list(x=s$v%*%(di*(t(s$u)%*%b)),q=sum(ii)))
}
