IORE.solution <- function(t, parent.0, k__iore, N)
{
	parent = (parent.0^(1 - N) - (1 - N) * k__iore * t)^(1/(1 - N))
}
