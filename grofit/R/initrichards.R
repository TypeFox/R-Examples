initrichards <-
function (time, y, A, mu, lambda)
{
if (is.numeric(time)==FALSE) stop("Need numeric vector for: time")
if (is.numeric(y)==FALSE) stop("Need numeric vector for: y")
if (is.numeric(mu)==FALSE) stop("Need numeric vector for: mu")
if (is.numeric(lambda)==FALSE) stop("Need numeric vector for: lambda")
if (is.numeric(A)==FALSE) stop("Need numeric vector for: A")

nu     <- 0.1
A      <- max(y)
mu     <- mu[1]
lambda <- lambda[1]

initrichards <- list (A=A, mu=mu, lambda=lambda, addpar=nu)
}

