#### functions for the F test in the univariate model
#### can be reused for selective inference too
#### slight degree of freedom adjustments if the intercept needs to be estimated
#### y is a matrix with the columns containing the replcates of the response
#### returns a list with:
####      - vector of F stat replications
####      - df1 and df2
compute.F.statistic <-
function(x, y, mu=NULL){
  ### precompute some quantities
  M = ncol(x)
  n = nrow(y)
  
  if (is.null(mu)){
    # X matrices
    X.tilde.1 = matrix (1, ncol=1, nrow=n)
    X.tilde.2 = cbind (X.tilde.1, x)
    
    # P matrix
    P.tilde.2 = X.tilde.2%*%solve(t(X.tilde.2)%*%X.tilde.2, t(X.tilde.2))
  }else{
    # P matrix
    P.tilde.2 = x%*%solve(t(x)%*%x, t(x))
  }
  
  ### compute the F stats
  
  df1 = M
  if (is.null(mu)){
    Py = P.tilde.2%*%y
    y1 = apply (y, 2, sum)
    yPy = apply(y*Py, 2, sum)
    yy = apply (y, 2, function(col){sum(col^2)})
    
    F.stats = (yPy - y1^2/n)*(n - M - 1)/(yy - yPy)/(M)
    
    # degrees of freedom
    df2 = n - M - 1
  }else{
    y.tilde = y - mu
    Py = P.tilde.2%*%y.tilde
    yPy = apply(y.tilde*Py, 2, sum)
    yy = apply (y.tilde, 2, function(col){sum(col^2)})
    
    F.stats = yPy*(n-M)/(yy - yPy)/M
    df2 = n - M
  }
  
  
  list(ts=F.stats, df1=df1, df2=df2)
}
