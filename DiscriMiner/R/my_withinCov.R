my_withinCov <-
function(X, g, div_by_n=FALSE)
{
  # X: matrix of explanatory variables
  # g: factor with group memberships
  # div_by_n: logical indicating division by num of observations

  # how many observations
  n = nrow(X)
  # how many variables
  ncx = ncol(X)
  # group levels and number of levels
  glevs = levels(g)
  ng = nlevels(g)
  # within cov matrix
  Within = matrix(0, ncx, ncx)
  for (k in 1:ng)
  {
    tmp <- g == glevs[k]
    nk = sum(tmp)
    if (div_by_n) {
      Wk = ((nk-1)/n) * var(X[tmp,])
    } else {
      #Wk = ((nk-1)/(n-1)) * var(X[tmp,])
      # divide by degrees of freedom
      Wk = ((nk-1)/(n-ng)) * var(X[tmp,])
    }
    Within = Within + Wk
    
  }
  # result
  Within
}
