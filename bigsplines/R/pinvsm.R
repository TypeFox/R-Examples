pinvsm <- 
function(X,iconst=-1,tol=.Machine$double.eps){
  # pseudoinverse for symmetric matrices
  
  xeig <- eigen(X,symmetric=TRUE)
  nze <- sum(xeig$val>xeig$val[1]*tol)
  Xinv <- xeig$vec[,1:nze]%*%((xeig$val[1:nze]^iconst)*t(xeig$vec[,1:nze]))
  
}