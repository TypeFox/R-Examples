my_discFunctions <-
function(X, g, group_means, within)
{
  # X: explanatory variables
  # g: factor with group memberships
  # group_means: group means matrix
  # within: pooled within-class covariance matrix
  
  # group means
  GM = group_means
  # how many groups
  ng = nrow(GM)
  # inverse of Within Cov Matrix
  W_inv = solve(within)
  
  # constant terms of fisher's discriminant linear functions
  alphas = rep(0, ng)
  # coefficients of fisher's discriminant linear functions
  betas = matrix(0, nrow(W_inv), ng)
  for (k in 1:ng) 
  {
    alphas[k] = -(1/2) * GM[k,] %*% W_inv %*% GM[k,]
    betas[,k] = t(GM[k,]) %*% W_inv
  }
  # Fisher's Discriminant Functions
  FDF = rbind(alphas, betas)
  rownames(FDF) = c("constant", colnames(X))
  colnames(FDF) = c(levels(g))
  # result
  FDF
}
