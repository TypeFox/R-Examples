dvolkov <- function(x, theta, m, J, order=96, log=FALSE){
	if(length(theta) > 1 | length(m) > 1 | length(J) > 1) stop ("Vectorization of parameters not implemented")
	if(!is.finite(theta) | theta <= 0 | !is.finite(J) | J <= 0 | m <= 0 | m >= 1) return (rep(NaN, length(x)))
	Cvolkov <- .C(volkov, res=as.double(rep(0,J)), theta0=as.double(theta), m0=as.double(m), J0=as.integer(J), N0=as.integer(order))
	y <- c(0, Cvolkov$res) # index 0 should receive 0
	if (any(!is.wholenumber(x))) warning("non integer values in x")
	x[ !is.wholenumber(x) | x < 0 | x > J ] <- 0
	y <- y[x+1]
	if (any(is.nan(y))) warning ("NaNs produced")
	if(log) return(log(y))
	else return(y)
}
