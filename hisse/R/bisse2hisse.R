#Simple function for converting output from diversitree to the proper HiSSE format in order to use reconstruction functions:

BisseToHisse <- function(lambda, mu, q01, q10){
	par.vec = rep(0, 56)
	par.vec[21:56] = 1
	par.vec[1:2] = lambda + mu
	par.vec[5:6] = mu / lambda
	par.vec[12] = q01
	par.vec[9] = q10
	return(par.vec)
}