my_quaDA <- 
function(X, y, learn, test, prior, prob)
{
  # Perform a quadratic discriminant analysis
  # X: matrix or data.frame with explanatory variables
  # y: vector or factor with group membership
  # learn: vector of learning observations
  # test: vector of testing observations
  # prior: vector of prior proportions
  # prob: logical indicating results in proability terms
  
  # how many observations
  n = nrow(X[learn,])
  ntest = length(test)
  # how many variables
  p = ncol(X)
  # how many groups
  ng = nlevels(y[learn])
  glevs = levels(y[learn])
  # how many obs in each group
  gobs = as.vector(table(y[learn]))
  names(gobs) = glevs
  if (any(gobs < p+1))
    stop("\nsome group category is too small for quaDA")
  
  # group means
  GM = my_groupMeans(X[learn,], y[learn])
  # within matrices based on QR decomposition
  WMqr = as.list(1:ng)
  # object to store log-determinants
  ldet = numeric(ng)
  # calculate ingredients
  for (k in 1:ng)
  {
    nk = gobs[k] - 1
    # center data
    Xcen = scale(X[unclass(y[learn])==k,], center=GM[k,], scale=FALSE)
    # QR decomposition
    qx = qr(Xcen / sqrt(nk) )
    if(qx$rank < p)
      stop("rank deficiency in group ", glevs[k])
    qx = qx$qr
    WMqr[[k]] = backsolve(qx[1:p, ], diag(p))
    ldet[k] = 2*sum(log(abs(diag(qx))))
  }
  
  ## classifcation
  # discrimination matrix to store results
  Disc <- matrix(0, nrow = ntest, ncol = ng)
  # calculate distances (the lower, the better)
  for (k in 1:ng) 
  {
    Xk = matrix(GM[k,], ntest, p, byrow = TRUE)
    dev = (X[test,] - Xk) %*% WMqr[[k]]
    Disc[,k] = 0.5 * rowSums(dev^2) + 0.5 * ldet[k] - log(prior[k])
  }
  
  # assignment in terms of probability?
  if (prob)
  {
    # Disc in terms of probability
    Disc <- exp( -(Disc - apply(Disc, 1, min, na.rm=TRUE))) 
    # predicting classes
    pred = Disc / drop(Disc %*% rep(1, ng))
    # predicted class
    pred_class = factor(max.col(pred), levels=seq_along(glevs), labels=glevs)
  } else {
    # predicted class
    pred = apply(Disc, 1, function(u) which(u == min(u)))
    names(pred) = NULL
    pred_class = factor(pred, levels=seq_along(glevs), labels=glevs)
  }
  dimnames(Disc) = list(rownames(X[test,]), glevs)
  # confusion matrix
  conf = table(original=y[test], predicted=pred_class)
  
  ## results
  res = list(WMqr = WMqr, GM = GM, ldet = ldet, prior = prior,
             Disc = Disc, 
             pred_class = pred_class, 
             conf = conf)
  res
}
