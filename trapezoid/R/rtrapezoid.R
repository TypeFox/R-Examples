rtrapezoid <- function(n, min = 0, mode1 = 1/3, mode2 = 2/3, max = 1, n1 = 2, n3 = 2, alpha = 1)
{
	out <- qtrapezoid(p = runif(n), min = min, mode1 = mode1, mode2 = mode2, max = max, n1 = n1, n3 = n3, alpha = alpha)
	return(out)
}
