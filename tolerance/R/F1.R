F1 <- function(a, b, b.prime, c, x, y, ...)
{
	A1.simple <- function(u, a, b, b.prime, c, x, y)
	{
		u^(a-1)*(1-u)^(c-a-1)*(1-u*x)^(-b)*(1-u*y)^(-b.prime)
	}
	gamma(c)/(gamma(a)*gamma(c-a))*as.numeric(integrate(A1.simple, 0, 1, a = a, b = b, b.prime = b.prime, c = c, x = x, y = y, ...)$value)
}

