my_mca <- 
function(X)
{
  # Perform multiple correspondence analysis
  # X: data frame with categorical variables as factors
  
  # how many observations
  nobs = nrow(X)
  # how many variables
  nvars = ncol(X)
  # number of categories per variable
  cats_per_var = sapply(X, nlevels)
  # total number of categories
  ncats = sum(cats_per_var)
  # number of factor coordinates (for MCA)
  nfacs = ncats - nvars
  
  # build super-indicator matrix Z
  Z = my_tdc(X)  
  # number of obs per category
  nopc = colSums(Z)
  # normalizing Z
  Znorm = sweep(Z, 2, sqrt(nvars*nopc), FUN="/")
  # apply svd
  Zsvd = svd(Znorm)
  # sequence with indices of components
  sec <- 1 + (1L:nfacs)
  # eigenvalues
  eigs = Zsvd$d[sec]^2
  values = cbind(eigs, 100*eigs/sum(eigs), 100*cumsum(eigs)/sum(eigs))
  colnames(values) = c("eigenvalues", "proportion", "accumulated")
  rownames(values) = 1:nfacs 
  # U-coefficients
  #U = diag(1/sqrt(nvars*nopc)) %*% Zsvd$v[,sec]/nvars
  U = diag(sqrt(nobs/(nvars*nopc))) %*% Zsvd$v[,sec]
  dimnames(U) = list(colnames(Z), paste("U",1:nfacs,sep=''))
  # row coordinates
  Fs = Z %*% U
  # add names
  dimnames(Fs) <- list(rownames(X), paste("F",1:nfacs,sep=''))  
  structure(list(values=values, coefficients=U, components=Fs), 
            class="qualmca")
}
