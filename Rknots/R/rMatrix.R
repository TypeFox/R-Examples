# Script comments and history
# 2011
# 5:35:25 PM

# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

rMatrix <-
function (theta) 
{
	if (missing(theta)) 
		stop("Argument theta missing, with no default\n")
	if (!is.numeric(theta))
	{
		warning("rotationMatrix: non numeric argument converted to numeric")
		theta <- as.numeric(theta)
	}
	R <- matrix(c(cos(theta), -sin(theta), 0, sin(theta), cos(theta), 
					0, 0, 0, 1), ncol = 3, byrow = TRUE)
	return(R)
}

