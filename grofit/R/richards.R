richards <-
function (time, A,mu,lambda,addpar)
{
A      <- A[1]
mu     <- mu[1]
lambda <- lambda[1]
nu     <- addpar[1]

if (is.numeric(time)==FALSE) stop("Need numeric vector for: time")
if (is.numeric(mu)==FALSE) stop("Need numeric vector for: mu")
if (is.numeric(lambda)==FALSE) stop("Need numeric vector for: lambda")
if (is.numeric(A)==FALSE) stop("Need numeric vector for: A")
if (is.numeric(nu)==FALSE) stop("Need numeric vector for: addpar[1]")

y        <- A*(1.0+nu*exp(1+nu)*exp(mu*(1+nu)^(1+1/nu)*(lambda-time)/A))^(-1/nu)
richards <- y
}

