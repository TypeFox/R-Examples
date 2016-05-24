# Inverts a matrix by using the QR decomposition and scaling
# Optionally imposes positive definiteness
.qrInverse <- function(mat, tolQR = 1e-7, imposePD = FALSE, imposeSym = FALSE, tilt = 0)
{
  
  if(imposePD){
    eigDec <- eigen( mat %*% mat + tilt * diag(nrow(mat)) )
    mat <- eigDec$vectors %*% ( t(eigDec$vectors) * sqrt(eigDec$value) )
  }
  
  D <- diag( diag(mat)^-0.5, nrow(mat) )
  
  inv <- D %*% qr.solve(D %*% mat %*% D, tol = tolQR) %*% D
  
  if( imposeSym ) inv <- ( inv + t(inv) ) / 2.0
  
  return( inv )
  
}