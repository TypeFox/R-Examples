resizeCircles <- function(x_only, y_only, overlap, standardDistance = sqrt(0.5))
{
	overlapArea <- function(r1, r2, d)
	{
		d1 <- (r1^2 - r2^2 + d^2) / (2 * d)
		d2 <- d - d1
		theta1 <- acos(d1 / r1)
		theta2 <- acos(d2 / r2)
		area1 <-  r1^2 * (theta1 - sin(2 * theta1) / 2)
		area2 <-  r2^2 * (theta2 - sin(2 * theta2) / 2)
		return (c(area1=area1, area2=area2, theta1=theta1, theta2=theta2, d1=d1, d2=d2))
	}

	if (x_only <= 0 | y_only <= 0 | overlap <= 0) return (NULL)
	
	ra1 <- sqrt((x_only + overlap) / pi)
	ra2 <- sqrt((y_only + overlap) / pi)
	
	r1 <- min(ra1,ra2)
	r2 <- max(ra1,ra2)
	bestFitDistance <- uniroot(function(d) sum(overlapArea(r1, r2, d)[1:2]) - overlap, c(r2 - r1+1e-6, r2 + r1 - 1e-6))$root

	return (c(ra1, ra2) / standardDistance / bestFitDistance)
}
