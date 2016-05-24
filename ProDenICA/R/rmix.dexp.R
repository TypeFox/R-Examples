rmix.dexp <-
function(n, means = 0)
{
	x <- logb(runif(n))
	x <- x * sample(c(-1, 1), n, replace = TRUE)
	xmeans <- if(length(means) > 1) sample(means, n, replace = TRUE) else 
			means
	x + xmeans
}

