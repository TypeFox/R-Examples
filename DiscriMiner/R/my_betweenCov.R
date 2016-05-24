my_betweenCov <-
function(X, g, div_by_n=FALSE)
{
  # X: matrix of explanatory variables
  # g: factor with group memberships
  
  # how many observations
  n = nrow(X)
  # how many variables
  p = ncol(X)
  # group levels and number of levels
  glevs = levels(g)
  ng = nlevels(g)
  # global mean
  mean_global = colMeans(X)
  # matrix to store results
  Between = matrix(0, p, p)
  # pooled between-class covariance matrix
  for (k in 1:ng)
  {
    # select obs of k-th group
    tmp <- g == glevs[k]
    # how many obs in group k
    nk = sum(tmp)
    # mean k-th group
    mean_k = colMeans(X[tmp,])
    # mean k-th group - global mean
    dif_k = mean_k - mean_global
    # k-th group between cov matrix
#    between_k = (nk/n) * tcrossprod(dif_k)
    if (div_by_n) {
      between_k = (nk/n) * tcrossprod(dif_k)
    } else {
      between_k = (nk/(n-1)) * tcrossprod(dif_k)
    }
    Between = Between + between_k
    
  }
  # result
  Between
}
