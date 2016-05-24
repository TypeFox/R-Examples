pbs <- function(q, N, S, lower.tail=TRUE, log.p = FALSE){
	N[ !is.finite(N) | N <= 0] <- NaN
	S[ !is.finite(S) | S <= 0] <- NaN
	y <- 1 + N*(1-q/N)^S/(q-N)
	y[q < 0] <- 0; y[q >= N] <- 1
	if(!lower.tail) y <- 1 - y
	if(log.p) y <- log(y)
	if(any(is.nan(y))) warning("NaNs produced")
	return(y)
}
