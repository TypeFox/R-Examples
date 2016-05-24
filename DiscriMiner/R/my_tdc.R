my_tdc <- 
function(X)
{  
  # Tableau Disjonctive Complete
  # (aka Complete Disjunctive Table)
  # X: data frame with categorical variables as factors
  
  # how many observations
  nobs = nrow(X)
  # how many variables
  nvars = ncol(X)
  # number of categories per variable
  cats_per_var = sapply(X, nlevels)
  # total number of categories
  ncats = sum(cats_per_var)
  
  # build super-indicator matrix Z
  Z = matrix(0, nobs, ncats)
  ini = cumsum(cats_per_var) - cats_per_var + 1
  fin = cumsum(cats_per_var)
  for (j in 1:nvars)
  {
    aux_lev = levels(X[,j])
    aux_mat = matrix(0, nobs, cats_per_var[j])
    for (k in 1:cats_per_var[j])
    {
      tmp <- X[,j] == aux_lev[k]
      aux_mat[tmp,k] = 1
    }
    Z[,ini[j]:fin[j]] = aux_mat
  }
  colnames(Z) = unlist(lapply(X, levels))
  Z
}
