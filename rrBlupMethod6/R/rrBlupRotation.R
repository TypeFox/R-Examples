rrBlupRotation <- function(y,
                           X = matrix(1,nrow = n, ncol = 1),
                           Z,
                           R)
{
  n <- length(y)
  
  svd_vare <- svd(R)
  
  W <- svd_vare$u %*% diag(svd_vare$d**-0.5) %*% t(svd_vare$v)
  
  y_tilda <- W %*% y
  Z_tilda <- W %*% Z
  X_tilda <- W %*% X
  
  return(list(y_tilda = y_tilda,
              X_tilda = X_tilda,
              Z_tilda = Z_tilda))
}
