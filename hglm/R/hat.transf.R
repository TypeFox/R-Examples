hat.transf <-
function(C22, transf, vc, w, k, N, tol.err = 1e-6) {
	qu <- numeric(k)
		middle <- diag(N) - C22/vc
		eigen.test <- eigen(middle, only.values = TRUE)
		#if (min(eigen.test$values) < tol.err) {
			## cat("Matrix middle bended\n")
			## cat("Minimum eigenvalue:\n", min(eigen.test$values))
		#	if (min(eigen.test$values) > (-tol.err)) middle <- middle + diag(N)*2*tol.err
		#	if (min(eigen.test$values) < (-tol.err)) middle <- middle - (min(eigen.test$values) - tol.err)*diag(N)
		#}
		M <- chol(middle)
		#if (!GPU) {
			qu <- 1 - colSums(tcrossprod(M, transf)**2)/w
		#}
		#else {
		#	qu <- 1 - colSums(gpuMatMult(M, t(transf))**2)/w
		#}
	return(qu)
}

