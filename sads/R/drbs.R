drbs <- function(x, N, S, log = FALSE) {
	N[ !is.finite(N) | N <= 0 | !is.wholenumber(N) ] <- NaN
	S[ !is.finite(S) | S <= 0 | !is.wholenumber(S) ] <- NaN
	y <- sapply(x, function(z)N/S*sum((z:S)^-1))/N
	if (any(is.nan(y))) warning ("NaNs produced")
	if (any(!is.wholenumber(x))) warning("non integer values in x")
	y[ ! is.wholenumber(x) | x < 1 | x > S ] <- 0
	if(!log) return(y)
	else return(log(y))
}
