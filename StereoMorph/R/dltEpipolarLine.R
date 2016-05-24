dltEpipolarLine <- function(p, cal.coeff1, cal.coeff2 = NULL, self=FALSE){
	# Translated and modified from Matlab function partialdlt written by Ty Hedrick
	# 'self' added, after Yekutieli et al. 2007

	# IF COEFFICIENT MATRIX IS INPUT TO FIRST COEFFICIENT PARAMETER
	if(is.null(cal.coeff2)){cal.coeff2 <- cal.coeff1[, 2];cal.coeff1 <- cal.coeff1[, 1];}

	# RANDOM Z VALUES (ACTUAL VALUES NOT IMPORTANT)
	z <- c(500,-500)

	# PREDICT X AND Y FOR ARBITRARY Z VALUES
	x <- z*NA
	y <- z*NA

	# PROJECT POINT P AT TWO Z VALUES
	for(i in 1:2){
		A <- p[1]*cal.coeff1[9]*cal.coeff1[7]*z[i] + p[1]*cal.coeff1[9]*cal.coeff1[8] - p[1]*cal.coeff1[11]*z[i]*cal.coeff1[5] -p[1]*cal.coeff1[5] + cal.coeff1[1]*p[2]*cal.coeff1[11]*z[i] + cal.coeff1[1]*p[2] - cal.coeff1[1]*cal.coeff1[7]*z[i] - cal.coeff1[1]*cal.coeff1[8] - cal.coeff1[3]*z[i]*p[2]*cal.coeff1[9] + cal.coeff1[3]*z[i]*cal.coeff1[5] - cal.coeff1[4]*p[2]*cal.coeff1[9] + cal.coeff1[4]*cal.coeff1[5]
		B <- p[1]*cal.coeff1[9]*cal.coeff1[6] - p[1]*cal.coeff1[10]*cal.coeff1[5] + cal.coeff1[1]*p[2]*cal.coeff1[10] - cal.coeff1[1]*cal.coeff1[6] - cal.coeff1[2]*p[2]*cal.coeff1[9] + cal.coeff1[2]*cal.coeff1[5]

		# SOLVE FOR X AND Y
		y[i] <- (-A) / B
		x[i] <- -(p[2]*cal.coeff1[10]*y[i] + p[2]*cal.coeff1[11]*z[i] + p[2]-cal.coeff1[6]*y[i] - cal.coeff1[7]*z[i] - cal.coeff1[8]) / (p[2]*cal.coeff1[9] - cal.coeff1[5])
	}

	# PROJECT 3D COORDINATE INTO SECOND IMAGE PLANE
	xy <- dltInverse(cal.coeff2, cbind(x,y,z))

	# GET THE LINE COEFFICIENTS
	m <- (xy[2,2] - xy[1,2]) / ((xy[2,1] - xy[1,1]) + 10^-10)
	b <- xy[1,2] - m*xy[1,1]

	# IF SPECIFIED, GET SELF-EPIPOLAR USING A POINT ON THE EPIPOLAR LINE
	if(self) return(dltEpipolarLine(p = c(0, b), cal.coeff2, cal.coeff1))

	# RETURN LINE COEFFICIENTS AND TWO POINTS ON THE LINE
	list(m = m, b = b, l1=c(0, b), l2=c(1, m + b))
}