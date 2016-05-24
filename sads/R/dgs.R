dgs <- function(x, k, S,  log = FALSE) {
	S[ !is.finite(S) | S < 1 | !is.wholenumber(S) ] <- NaN
	k[ k<0 | k>1 ] <- NaN
	cf <- 1/(1-(1-k)^S)
	y <- cf*k*(1-k)^(x-1)
	if (any(is.nan(y))) warning ("NaNs produced")
	if (any(!is.wholenumber(x))) warning("non integer values in x")
	y[ ! is.wholenumber(x) | x < 1] <- 0
	if(!log) return(y)
	else return(log(y))
}
