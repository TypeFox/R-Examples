my_geoDA <- 
function(X, y, learn, test)
{
  # Perform a geometric predictive discriminant analysis
  # X: matrix or data.frame with explanatory variables
  # y: vector or factor with group membership
  # learn: vector of learning observations
  # test: vector of testing observations
  
  # how many observations
  n = nrow(X[learn,])
  ntest = length(test)
  # how many groups
  ng = nlevels(y[learn])
  glevs = levels(y[learn])
  # group means
  GM = my_groupMeans(X[learn,], y[learn])
  # within-class covariance matrix
  W = my_withinCov(X[learn,], y[learn])
  # inverse of Within cov matrix
  W_inv = solve(W)
  
  # constant terms of classification functions
  alphas = rep(0, ng)
  # coefficients of classification functions
  Betas = matrix(0, nrow(W_inv), ng)
  for (k in 1:ng) 
  {
    alphas[k] = -(1/2) * GM[k,] %*% W_inv %*% GM[k,]
    Betas[,k] = t(GM[k,]) %*% W_inv
  }
  # Mahalanobis-Fisher Classification Rule
  FDF = rbind(alphas, Betas)
  rownames(FDF) = c("constant", colnames(X))
  colnames(FDF) = glevs 
  
  # matrix of constant terms
  A = matrix(rep(alphas,ntest), ntest, ng, byrow=TRUE)
  # apply discrim functions
  Disc = X[test,] %*% Betas + A
  dimnames(Disc) = list(rownames(X[test,]), glevs)
  # predicted class
  pred = apply(Disc, 1, function(u) which(u == max(u)))
  names(pred) = NULL
  # assign class values
  pred_class = factor(pred, levels=seq_along(glevs), labels=glevs)
  
  # confusion matrix
  conf = table(original=y[test], predicted=pred_class)
  # results
  res = list(FDF=FDF, conf=conf, Disc=Disc, pred_class=pred_class)
}
