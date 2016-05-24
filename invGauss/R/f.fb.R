"f.fb" <-
function(t, mu = 0.341, tau = 0.247, c1 = 2.58)
{
# COMPUTES (DEFECTIVE) DENSITY, FORMULA (7)
#
	c1/(sqrt(2 * pi) * t * sqrt(t^2 * tau^2 + t)) * exp( - (c1 - mu * t)^2/(
		2 * (t^2 * tau^2 + t)))
}

