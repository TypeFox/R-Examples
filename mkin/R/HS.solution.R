HS.solution <- function(t, parent.0, k1, k2, tb)
{
	parent = ifelse(t <= tb, 
		parent.0 * exp(-k1 * t),
		parent.0 * exp(-k1 * tb) * exp(-k2 * (t - tb)))
}
