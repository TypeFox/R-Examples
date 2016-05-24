dmzsm <- function(x, J, theta, log = FALSE){
	J[ !is.finite(J) | J <= 0] <- NaN
	theta[ !is.finite(theta) | theta <= 0] <- NaN
	mzsm <- function(y, J, theta) (theta/y)*(1-y/J)^(theta-1)
	sn <- mzsm(y=x, J = J, theta = theta)
	mu <- mzsm(y=1:J, J = J, theta = theta)
	lpn <- suppressWarnings(log(sn) - log(sum(mu))) # Might be NaN if x > 0, but end result will be zero
        lpn[ x <= 0 | x > J ] <- -Inf
	if (any(is.nan(lpn))) warning ("NaNs produced")
	if(log) return(lpn)
	else return(exp(lpn))
}
