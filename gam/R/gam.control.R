"gam.control" <-
function(epsilon = 9.9999999999999995e-08, bf.epsilon = 9.9999999999999995e-08,
	maxit = 30, bf.maxit = 30, trace = FALSE, ...)
{
	if(epsilon <= 0) {
		warning("the value of epsilon supplied is zero or negative; the default value of 1e-7 was used instead"
			)
		epsilon <- 9.9999999999999995e-08
	}
	if(maxit < 1) {
		warning("the value of maxit supplied is too small; the default value of 30 was used instead"
			)
		maxit <- 30
	}
	if(bf.epsilon <= 0) {
		warning("the value of bf.epsilon supplied is zero or negative; the default value of 1e-7 was used instead"
			)
		bf.epsilon <- 9.9999999999999995e-08
	}
	if(bf.maxit < 1) {
		warning("the value of bf.maxit supplied is too small; the default value of 30 was used instead"
			)
		bf.maxit <- 30
	}
	list(epsilon = epsilon, maxit = maxit, bf.epsilon = bf.epsilon, 
		bf.maxit = bf.maxit, trace = as.logical(trace)[1])
}
