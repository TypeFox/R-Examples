msBP.leaf.allocation <-
function(y, scale=4)
{
n <- length(y)
s <- rep(scale, n)
h <- rep(NA, n)
for(i in 1:n)
{
	prob <- rep(0, 2^scale)
	for(j in 1:(2^scale))
	{
		prob[j] <- dbeta(y[i], j, 2^scale - j +1)
	}
	h[i] <- min(c(1:(2^scale))[prob == max(prob)])
}
list(s=s, h=h)
}
