dmand <- function (x, N, s, v, log = FALSE) {
	N[ !is.finite(N) | N < 1 | !is.wholenumber(N) ] <- NaN
	s[ !is.finite(s) | s <= 0] <- NaN
	v[ !is.finite(v) | v < 0] <- NaN
	lny <- - s * log(x+v) - log(sum(((1:N)+v)^(-s)))
	if (any(is.nan(lny))) warning ("NaNs produced")
	if (any(!is.wholenumber(x))) warning("non integer values in x")
	lny[ ! is.wholenumber(x) | x < 1] <- -Inf
	if (log) return(lny)
	else return(exp(lny))
}
