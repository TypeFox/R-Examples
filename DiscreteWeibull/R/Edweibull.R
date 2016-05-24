Edweibull<-function (q, beta, eps = 1e-04, nmax = 1000, zero = FALSE) 
{
    if (beta == 1 & !zero) 
        e <- 1/(1 - q)

    else if (beta == 1 & zero) 
        e <- q/(1 - q)
    else
	{
        xmax <- min(2 * qdweibull(1 - eps, q, beta, zero), nmax)
	if(xmax < nmax)
	{
        x <- 1:xmax
        e<-sum(ddweibull(x, q, beta, zero) * x)
	}
	else # approximation with expected value of continuous Weibull
	{
	lambda<-(-1/log(q))^(1/beta)
	e <- lambda*gamma(1+1/beta)
	e <- e + 1/2 -zero 
        }
	}
	return(e)
}


E2dweibull<-function (q, beta, eps = 1e-04, nmax = 1000, zero = FALSE) 
{
    if (beta == 1 & !zero) 
        e <- (1 + q)/(1 - q)^2
    if (beta == 1 & zero) 
        e <- q * (1 + q)/(1 - q)^2
    else {
        xmax <- 2*qdweibull(1 - eps, q, beta, zero)
        if (xmax < nmax)
	{
        x <- 1:xmax
        e <- sum(ddweibull(x, q, beta, zero) * x^2)
	}
	else # approximation
	{
	lambda <- (-1/log(q))^(1/beta)
	e <- lambda^2*(gamma(1+2/beta)-(gamma(1+1/beta)^2)) +
	(Edweibull(q, beta, eps = eps, nmax = nmax, zero = zero))^2
	e <- ceiling(e) - zero
	}
        }
    return(e)
}

Vdweibull <-
function(q, beta, eps=0.0001, nmax=1000, zero=FALSE)
{
if(beta==1)
q/(1-q)^2
else
{
E2dweibull(q, beta, eps, nmax, zero) - Edweibull(q, beta, eps, nmax, zero)^2
}
}

ERdweibull <-
function(q, beta, eps=0.0001, nmax=1000)
{
if(beta==1)
(1-q)/q*log(1/(1-q))
else
{
xmax<-min(2*qdweibull(1-eps, q, beta),nmax)
x <- 1:xmax
sum((q^(x-1)^beta-q^x^beta)*1/x)
}
}

