fpca.start <- function(A, maxsteps = 200, tol = 1e-3, normalised = T, K = 2, score = F, ridge = T, approx = F){

	if(score == F) return.obj = fpca.nonscore(A = A, maxsteps = maxsteps, tol = tol, normalised = normalised, K = K, ridge = ridge, approx = approx)
	
	if(score == T) return.obj = fpca.score(A = A, maxsteps = maxsteps, tol = tol, K = K, ridge = ridge, approx = approx)
	
	return(return.obj)
}
