`sincdata` <-
function(n, noise = 0, a = -4, b = 4)
{
	# test if n is a number and > 0
	if(length(n) != 1 || n <= 0)
	{
		print("n > 0 must be the number of samples to generate")
	}
	# test if noise is a number
	if(length(noise) != 1)
	{
		print("noise must be a number specifying the noise level")
	}
	# test if 'a' and 'b' are numbers and a < b
	if(length(a) != 1 || length(b) != 1 || a >= b)
	{
		print("a and b must be real numbers with a < b specifying 
		the interval [a,b] from which ths xs are drawn")
	}
	
	# draw n points uniformly from [a, b]
	X <- matrix(runif(n, a, b), n, 1)
	
	# calculate y's + normal-noise
	y <- matrix(sinc(X) + noise*rnorm(n), 1)
	
	return(list(X = X, y = y))
}

