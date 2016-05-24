rarefy <-
function(raw, fraction)
{
	if (fraction > 1)
	{
		stop("fraction must be less than 1")	
	}
	
	nr <- nrow(raw)
	
	indices <- sample(nr, floor(fraction*nr))
	
	return(raw[indices,])	
}

