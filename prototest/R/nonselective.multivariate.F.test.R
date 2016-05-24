#### computes the multivarite F statistic
#### testing the omission of the 'test.group' columns labelled by 'groups'
#### y is a MATRIX of replications of the response vector (each column a different replication)
nonselective.multivariate.F.test <-
function(x, y, groups, test.group, mu){
  ## extract some quantities
  n = nrow (y)
  p = ncol(x)
  M = sum (groups == test.group)
  
  ## do we estimate the intercept?
  X1 = X2 = NULL
  if (is.null(mu)){
    X1 = cbind (X1, rep(1, n))
    X2 = cbind (X2, rep(1, n))
  }
  X2 = cbind (X2, x)
  P2 = X2%*%ginv(X2)
  df2 = sum(diag(P2))
  
  ## response and mu
  y.tilde = y
  if (!is.null(mu)){
    y.tilde = y.tilde-mu
  }
  P2y = P2%*%y.tilde
  yP2y = apply(y.tilde*P2y, 2, sum)
  yy = apply (y.tilde^2, 2, sum)
  
  ## find the two projection matrices
  non.test.selector = groups != test.group
  if (sum(non.test.selector) == 0){ # no columns from the non-interest groups
    if (is.null(X1)){# no mean estimated, so no yP1y in numerator
      F.stat = yP2y/(yy - yP2y)/df2*(n - df2)
      df1 = 0
    }else{
      ## other projection matrix in numerator
      P1 = X1%*%ginv(X1)
      df1 = sum(diag(P1))
      P1y = P1%*%y
      yP1y = apply(y*P1y, 2, sum)
      
      F.stat = (yP2y - yP1y)/(yy - yP2y)/(df2-df1)*(n-df2)
    }
  }else{ ## columns selected from other 
    ## projection matrix
    X1 = cbind(X1, x[,non.test.selector, drop=FALSE])
    P1 = X1%*%ginv(X1)
    df1 = sum(diag(P1))
    
    ## multiply with y
    P1y = P1%*%y
    yP1y = apply(y*P1y, 2, sum)
    
    F.stat = (yP2y - yP1y)/(yy - yP2y)/(df2-df1)*(n-df2)
  }
  
  p.val = pf (F.stat, df1=df2-df1, df2=n-df2, lower.tail=FALSE)
  
  return (list(ts=F.stat, p.val=p.val))
}
