cronbach.alpha <- function (X)
{
	#---------------------------------------------------
	# Cronbach's alpha (standardized version)
	#---------------------------------------------------
	nvar <- ncol(X)
	crX <- cor(X)
	sumcorr <- sum(sum(crX[upper.tri(crX, diag = FALSE)]))
	alpha <- nvar*sumcorr/(.5*nvar*(nvar-1)+(nvar-1)*sumcorr)
	
	alpha
}