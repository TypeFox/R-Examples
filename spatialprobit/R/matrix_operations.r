require(Matrix)

# PURPOSE: Creates an nxn sparse identity matrix 
#---------------------------------------------------
# USAGE: result = speye(n)
#---------------------------------------------------
# RETURNS: a sparse n x n identity matrix
# --------------------------------------------------
speye <- function(n){
  In <- Matrix(data=0,ncol=n,nrow=n,sparse=T)
  diag(In) <- 1
  return(In)
}

# PURPOSE: converts a three column matrix in a 
#          sparse matrix
#---------------------------------------------------
# USAGE: result = spconvert(X)
# where: X is a three colum matrix that should be 
#        interpreted as [ nrow ncol value ]
#---------------------------------------------------
# RETURNS: A  martix where W[nrow,ncol] = value
#          as defined by each row of the X matrix
# --------------------------------------------------
spconvert <- function(X){
  n      <- max(X[,1])
  m      <- max(X[,2])
  Xrlist <- which( X[,3] != 0 )
  result <- Matrix(data=0,nrow=n,ncol=m,sparse=T)
  for(Xr in Xrlist ){
    result[ X[Xr,1], X[Xr,2] ] <- X[Xr, 3]
  }
  return(result)
}

