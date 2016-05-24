my_linDA <- 
function(X, y, learn, test, prior, prob)
{
  # Perform a lienar discriminant analysis
  # X: matrix or data.frame with explanatory variables
  # y: vector or factor with group membership
  # learn: vector of learning observations
  # test: vector of testing observations
  # prior: vector of prior proportions
  # prob: logical indicating results in proability terms
  
  # how many observations
  n = nrow(X[learn,])
  ntest = length(test)
  # how many groups
  ng = nlevels(y[learn])
  glevs = levels(y[learn])
  # how many obs in each group
  nobs_group = as.vector(table(y[learn]))
  # group means
  GM = my_groupMeans(X[learn,], y[learn])
  # within-class covariance matrix
  W = my_withinCov(X[learn,], y[learn])
  # inverse of Within cov matrix
  W_inv = solve(W)
  
  # constant terms of classification functions
  cons = rep(0, ng)
  # coefficients of classification functions
  Betas = matrix(0, nrow(W_inv), ng)
  for (k in 1:ng) 
  {
    cons[k] = -(1/2) * GM[k,] %*% W_inv %*% GM[k,] + log(prior[k])
    Betas[,k] = t(GM[k,]) %*% W_inv
  }
  # Fisher's Discriminant Functions
  FDF = rbind(cons, Betas)
  rownames(FDF) = c("constant", colnames(X))
  colnames(FDF) = glevs 
  
  # matrix of constant terms
  A = matrix(rep(cons,ntest), ntest, ng, byrow=TRUE)
  # apply discrim functions
  Disc = X[test,] %*% Betas + A
  
  # probability values 
  if (prob) {
    # exponential
    Disc <- 1 - exp( -(Disc - apply(Disc, 1, min, na.rm=TRUE)))
    # predicting classes
    pred = Disc / drop(Disc %*% rep(1, ng))
    # predicted class
    pred_class = factor(max.col(pred), levels=seq_along(glevs), labels=glevs)
  } else {
    # predicted class
    pred = apply(Disc, 1, function(u) which(u == max(u)))
    names(pred) = NULL
    # assign class values
    pred_class = factor(pred, levels=seq_along(glevs), labels=glevs)    
  }
  dimnames(Disc) = list(rownames(X[test,]), glevs)
  # confusion matrix
  conf = table(original=y[test], predicted=pred_class)
  # results
  res = list(FDF=FDF, conf=conf, Disc=Disc, pred_class=pred_class)
  res
}
