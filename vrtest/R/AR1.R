AR1 <-
function(x)
{
T <- length(x) -1
Y <- x[2:T]
Y_ <- x[1:(T-1)]

ALPHA <- solve(t(Y_) %*% Y_ ) %*% t(Y_) %*% Y
RE <- Y - Y_ %*% ALPHA
SIGMAS <- sum(RE^2) / (T-1)
STDA <- sqrt( SIGMAS * solve(t(Y_) %*% Y_ ))
return(list(ALPHA=ALPHA,STDA=STDA))
}
