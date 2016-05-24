initgompertz.exp <-
function (time, y, A, mu, lambda)
{

if (is.numeric(time)==FALSE) stop("Need numeric vector for: time")
if (is.numeric(y)==FALSE) stop("Need numeric vector for: y")
if (is.numeric(mu)==FALSE) stop("Need numeric vector for: mu")
if (is.numeric(lambda)==FALSE) stop("Need numeric vector for: lambda")
if (is.numeric(A)==FALSE) stop("Need numeric vector for: A")

alfa   <- 0.1
tshift <- max(time)/10
A      <- max(y)
mu     <- mu[1]
lambda <- lambda[1]

initgompertz.exp <- list (A=A, mu=mu, lambda=lambda, addpar=c(alfa, tshift))
}

