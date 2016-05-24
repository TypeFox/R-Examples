`logspace` <-
function(l, u, n)
{
	# test if u and l are numbers and l <= u
	if(length(l) != 1 || length(u) != 1 || l > u)
	{
		print("l and u must be numbers and l <= u")
		return()
	}
	# test if na is a number
	if(length(n) != 1)
	{
		print("n must be the number of points to generate")
		return()
	}
	
	return(10^seq(l, u, length.out=n))
}

