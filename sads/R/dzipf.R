dzipf <- function(x, N, s, log = FALSE) {
	N[ !is.finite(N) | N <= 0 | !is.wholenumber(N) ] <- NaN
	s[ !is.finite(s) | s <= 0] <- NaN
	y <- -s*log(x)-log(sum(1/(1:N)^s))
	if (any(is.nan(y))) warning ("NaNs produced")
	if (any(!is.wholenumber(x))) warning("non integer values in x")
	y[ ! is.wholenumber(x) | x < 1 | x > N] <- 0
	if(log) return(y)
	else return(exp(y))
}
