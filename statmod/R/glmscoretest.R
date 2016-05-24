##  glmscore.R

glm.scoretest <- function(fit, x2, dispersion=NULL)
#	Score test for new covariate in glm
#	Gordon Smyth
#	27 March 2009. Last modified 20 Mar 2010.
{
	w <- fit$weights
	r <- fit$residuals
	if(any(w <= 0)) {
		r <- r[w>0]
		x2 <- x2[w>0]
		w <- w[w>0]
	}
	if (is.null(dispersion)) {
		fixed.dispersion <- (fit$family$family %in% c("poisson","binomial")) 
		if(fixed.dispersion)
			dispersion <- 1
		else if(fit$df.residual > 0) {
			dispersion <- sum(w*r^2)/fit$df.residual
		} else {
			stop("No residual df available to estimate dispersion")
		}
	}
	ws <- sqrt(w)
	x2.1w <- qr.resid(fit$qr,ws*x2)
	zw <- ws*r
	colSums(as.matrix(x2.1w*zw))/sqrt(colSums(as.matrix(x2.1w * x2.1w)))/sqrt(dispersion)
}

