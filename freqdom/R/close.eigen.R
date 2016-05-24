# NOT EXPORTED
# While computing eigenspaces in complex setup we deal with a problem that
# eigendirections are not unique - they can be rotated by any scalar \eqn{z} in \eqn{C}, s.t. \eqn{|z| = 1}
# In this case we have an inifinite number of versions. Function \code{close.eigen}
# finds the eigenvalues of \code{M} minimizing the \eqn{L_2} distance from \code{Prev}.
#
# @title Find close eigendirections
# @param M matrix to compute close eigenvalues from rotate
# @param Prev matrix to approximate
# @return rotated matrix \code{M}
# @importFrom stats rnorm
# @export
# @examples
# M = matrix(rnorm(9),3,3) + matrix(rnorm(9),3,3)*1i
# M = M%*%t(M)
# E1 = eigen(M + diag(3) * 0.1)
# E = close.eigen(M,E1$vectors)
close.eigen = function(M,Prev){
  if (!is.matrix(M) || dim(M)[1] != dim(M)[2])
    stop ("M must be a square matrix")
  if (!is.matrix(Prev) || dim(Prev)[1] != dim(Prev)[2])
    stop ("Prev must be a square matrix")
  if (dim(Prev)[1] != dim(M)[2])
    stop ("Dimensions of M and P must be equal")
  
  Eg = eigen(M)
	V = Re(Eg$vectors)
	W = Re(Prev)
	nbasis = dim(M)[1]

	for (col in 1:nbasis){		
		if (sum(V[,col] * W[,col]) < 0)
			Eg$vectors[,col] = -Eg$vectors[,col]
	}

	Eg
}
