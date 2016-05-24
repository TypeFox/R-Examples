hatfun <- function(x,g)
{
	Q  <- x%*%solve(t(x)%*%(g*x))%*%t(g*x)
	1 + diag(t(Q)%*%Q - Q - t(Q))
}
