my_catDA <- 
function(X, y, learn, test, autosel, prob)
{
  # Perform a geometric predictive discriminant analysis
  # X: data frame with categorical explanatory variables
  # y: vector or factor with group membership
  # learn: vector of learning observations
  # test: vector of testing observations
  # autosel: logical indicating automatic selection of MCA comps
  # prob: probability level for automatic selection
  
  ## main ingredients
  # how many observations
  n = nrow(X[learn,])
  ntest = length(test)
  # how many variables
  p = ncol(X)
  # how many groups
  ng = nlevels(y[learn])
  glevs = levels(y[learn])
  # MCA
  getmca = my_mca(X[learn,])
  # number of factors
  nf = min(p, ng-1)
  # selection of components
  select = 1:ncol(getmca$components)
  if (autosel) {
    getpower = discPower(getmca$components, y[learn])
    select = which(getpower$p_value <= prob)
  }
  # coordinate factors (MCA components)
  Fs = getmca$components[,select]
  # mca coefficients
  U = getmca$coefficients[,select]
  # super-indicator matrix (TDC)
  Z = my_tdc(X[learn,])
  
  ## geometric discriminant analysis
  # group means
  GM = my_groupMeans(Fs, y[learn])
  # within-class covariance matrix
  W = my_withinCov(Fs, y[learn])
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
  dimnames(Betas) = list(colnames(getmca$components[,select]), glevs)
  # raw coeffs of classification functions
  Raw = U %*% Betas
  
  # transforming coefficients 0-1000
  cats_per_var = sapply(X[learn,], nlevels)
  fin = cumsum(cats_per_var)
  ini = fin - cats_per_var + 1
  # normalized coefficients
  Norm = matrix(0, nrow(Raw), ncol(Raw))
  for (k in 1:ng)
  {
    for (j in 1:p)
    {
      tmp = Raw[ini[j]:fin[j],k]
      Norm[ini[j]:fin[j],k] = tmp - min(tmp)
    }
    Norm[,k] = Norm[,k] * (1000/sum(Norm[,k]))
  }
  dimnames(Norm) = dimnames(Raw)
  
  ## classification
  # (there is no constant term, Saporta 2006, page 462)
  # apply discrim functions
  Ztest = my_tdc(X[test,])
  Disc = Ztest %*% Norm
  # predicted class
  pred = apply(Disc, 1, function(x) which(x == max(x)))
  # assign class values
  pred_class = factor(pred, levels=seq_along(glevs), labels=glevs)
  
  # confusion matrix
  conf = table(original=y[test], predicted=pred_class)
  # results
  res = list(Raw = Raw, 
             Norm = Norm, 
             conf = conf, 
             Disc = Disc, 
             pred_class = pred_class)
}
