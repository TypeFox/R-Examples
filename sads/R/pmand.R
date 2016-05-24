pmand <- function (q, N, s, v, lower.tail = TRUE, log.p = FALSE){
	N[ !is.finite(N) | N < 1 | !is.wholenumber(N) ] <- NaN
	s[ !is.finite(s) | s <= 0] <- NaN
	v[ !is.finite(v) | v < 0] <- NaN
	y <- c()
	for (i in 1:length(q)) {
		if (is.nan(q[i])) y[i] <- NaN
		else y[i] <- log(sum(1/((1:q[i])+v)^s)) - log(sum(1/((1:N)+v)^s))
	}
    y <- exp(y)
	if (any(!is.wholenumber(q))) warning("non integer values in q")
	y[ ! is.wholenumber(q) | q < 1 ] <- 0
    if (!lower.tail) y <- 1 - y
    if (log.p) y <- log(y)
	if (any(is.nan(y))) warning ("NaNs produced")
  return(y)
}
