sqrtmatrix <-
function(mat)
{
  # Checking the symmetry of the matrix
	if (!isSymmetric(mat))
    stop("The matrix must be symmetric")
  
  ep<-eigen(mat,symmetric=TRUE);
  # Checking the positivity of the eigenvalues
  i.neg=which(ep$values < -.Machine$double.eps)
	if (length(i.neg) > 0)
    {stop("All eigenvalues of the matrix must be greater or equal to zero")}
  
  # Transforming the small eigenvalues to zero
  i.small<-which(abs(ep$values) < .Machine$double.eps)
  if (length(i.small)>0)
     {ep$values[i.small] = 0}

	valp = diag(sqrt(abs(ep$values)));
	vecp = ep$vectors;
	vecp %*% valp %*% t(vecp)
}
