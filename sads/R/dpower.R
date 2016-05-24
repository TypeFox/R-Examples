dpower <- function(x, s, log = FALSE){
	# Circumventing errors in function zeta
	to.NaN <- !is.finite(s) | s <= 1
	s[to.NaN] <- 0
	y <- -s*log(x)-log(zeta(s))
	y[to.NaN] <- NaN
	if (any(is.nan(y))) warning ("NaNs produced")
	if (any(!is.wholenumber(x))) warning("non integer values in x")
	y[ ! is.wholenumber(x) | x < 1] <- -Inf
	if(log) return(y)
	else return(exp(y))
}
