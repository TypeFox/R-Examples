dtrapezoid <- function(x, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, 
		n3 = 2, alpha = 1, log = FALSE)
{
	# Calculate density from C function
	out <- .Call("dtrapezoid", x, min, mode1, mode2, max, n1, n3, alpha)
	
	if(log) out <- log(out)
	if(any(is.nan(out))) 
	{
		warning("NaN in dtrapezoid")
	} else if(any(is.na(out))) 
	{
		warning("NA in dtrapezoid")
	}
	
	return(out)
}
