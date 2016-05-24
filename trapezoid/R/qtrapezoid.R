qtrapezoid <- function(p, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, 
		n3 = 2, alpha = 1, lower.tail = TRUE, log.p = FALSE)
{
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1-p
	
	# Calculate quantiles from C function
	out <- .Call("qtrapezoid", p, min, mode1, mode2, max, n1, n3, alpha)
	
	return(out)
}
