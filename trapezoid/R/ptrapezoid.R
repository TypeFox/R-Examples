ptrapezoid <- function(q, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, 
		n3 = 2, alpha = 1, lower.tail = TRUE, log.p = FALSE)
{
	# Calculate cumulative probability from C function
	out <- .Call("ptrapezoid", q, min, mode1, mode2, max, n1, n3, alpha)
	
	if(!lower.tail) out <- 1 - out
	if(log.p) out <- log(out)
	if(any(is.nan(out))) 
	{
		warning("NaN in dtrapezoid")
	} else if(any(is.na(out))) 
	{
		warning("NA in dtrapezoid")
	}
	
	return(out)	
}
