gompertz <-
function (time, A,mu,lambda, addpar=NULL)
{
A      <- A[1]
mu     <- mu[1]
lambda <- lambda[1]

if (is.numeric(time)==FALSE) stop("Need numeric vector for: time")
if (is.numeric(mu)==FALSE) stop("Need numeric vector for: mu")
if (is.numeric(lambda)==FALSE) stop("Need numeric vector for: lambda")
if (is.numeric(A)==FALSE) stop("Need numeric vector for: A")

e        <- exp(1)
y        <- A*exp(-exp(mu*e*(lambda-time)/A +1.0))
gompertz <- y
}

