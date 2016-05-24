DFOP.solution <- function(t, parent.0, k1, k2, g)
{
	parent = g * parent.0 * exp(-k1 * t) +
		 (1 - g) * parent.0 * exp(-k2 * t)
}
