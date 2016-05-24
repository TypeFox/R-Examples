`calcLogLikeManual` <-
function(gp) {
	N = dim(as.matrix(gp$Z))[1]
	ans = -N / 2 * log(2*pi) - .5 * determinant(solve(gp$invVarMatrix), logarithm=TRUE)$modulus[1] 
	ans = ans - .5*t(gp$Z - gp$mu)%*%gp$invVarMatrix%*%(gp$Z - gp$mu)
	return (ans)
}

