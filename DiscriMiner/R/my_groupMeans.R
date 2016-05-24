my_groupMeans <-
function(X, g)
{
  # X: matrix of explanatory variables
  # g: factor with group memberships
  
  # how many groups
  ng = nlevels(g)
  # matrix with group means
  Means = matrix(0, ng, ncol(X))
  for (j in 1:ncol(X))
  {
    Means[,j] = tapply(X[,j], g, FUN=mean)
  }
  # add names
  rownames(Means) = levels(g)
  colnames(Means) = colnames(X)
  # results
  Means
}
