com.log.sum = function(x,y)		# log.sum(x,y) = log( exp(x) + exp(y) )
{
	if (x == -Inf)
		{ return (y); }
	else if (y == -Inf)
		{ return (x); }
	else if (x > y)
		{ return (x + log( 1 + exp(y - x) ) ); }
	else
		{ return (y + log( 1 + exp(x - y) ) ); }
}

com.log.difference = function(x,y)	# log.difference(x,y) = log( exp(x) - exp(y) )
{
	if (x == -Inf)
		{ return (NaN); }
	else if (y == -Inf)
		{ return (x); }
	else if (x > y)
		{ return (x + log( 1 - exp(y - x) ) ); }
	else
		{ return (NaN); }
}


com.log.factorial = function(x)	# log(factorial(x))
{
	if (is.vector(x) && length(x) > 1)
	{
		for (i in 1:length(x))
			x[i] = com.log.factorial(x[i]);
		return (x);
	}
	else if (is.numeric(x))
	{
		if (x == 0) { x = 1; }
		return (sum(log(seq(from = 1, to = x, by = 1))));
	}
	else { stop("x must be a vector or number."); }
}
