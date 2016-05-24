FOMC.solution <- function(t, parent.0, alpha, beta)
{
	parent = parent.0 / (t/beta + 1)^alpha
}
