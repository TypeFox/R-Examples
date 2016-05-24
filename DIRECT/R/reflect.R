reflect <-
function (x, lower, upper)
{
	L = upper - lower
	if (x>=lower & x<=upper)
		result = x
	else
	{
		if (x<lower)
		{
			near = lower
			far = upper
		}
		if (x>upper)
		{
			near = upper
			far = lower
		}
	
		# compute distance between x and nearest boundary
		d = abs (near - x)
		# compute number of intervals contained in this distance
		n = floor (d/L)
		# n odd: subtract remainder from upper boundary
		# n even: add remainder to lower boundary
		d.tmp = d - n*L
		if (x<lower)
			result = ifelse ((n%%2)!=0, far-d.tmp, near+d.tmp)
		else
			result = ifelse ((n%%2)!=0, far+d.tmp, near-d.tmp)
	}
	
	return (result)
}

