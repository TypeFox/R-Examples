.savGol <-
function(ndvi, window.sav=7, degree=2, smoothing=10){
	if (!is.numeric(window.sav)){ stop("'window.sav' should be of type 'numeric'") }	
	if (!is.numeric(degree)){ stop("'degree' should be of type 'numeric'") }
	if (!is.numeric(smoothing)){ stop("'smoothing' should be of type 'numeric'") }
	createCoef <- function(window, degree){
		degree <- degree + 1
		A <- matrix(NA, window, degree)
		base <- (1:window - ceiling(window/2))
		exp <- 0:(degree-1)
		for (i in 1:window){
			for (j in 1:degree){
				A[i,j] <- base[i]^exp[j]
			}
		}
		Inv <- qr.solve(t(A) %*% A)
		E <- matrix(0, window, window)
		for (i in 1:window){
			E[i,i] <- 1
		}
		c <- vector(mode="numeric", length=window)
		for (i in 1:window){
			c[i] <- (Inv %*% (t(A) %*% E[i,]))[1]
		}
		return(c)
	}

	#initialize
	days <- length(ndvi)
	ndvi <- ifelse(is.na(ndvi), -1, ndvi)
	correctedndvi <- vector(mode="numeric", length=days)

	#Savitzky-Golay Smoothing
	coef <- vector(mode="numeric", length=window.sav);
	coef <- createCoef(window.sav, degree);

	res <- .C("savGol", rdays=as.integer(days), ndvi=as.numeric(ndvi), 
		coef=as.numeric(coef),nCoef=as.integer(window.sav) ,cndvi=as.numeric(correctedndvi),
		PACKAGE="phenex")$cndvi

	for (i in 0:(smoothing-1)){	
		res <- .C("savGol", rdays=as.integer(days), ndvi=as.numeric(res), 
			coef=as.numeric(coef),nCoef=as.integer(window.sav) ,cndvi=as.numeric(correctedndvi),
			PACKAGE="phenex")$cndvi
	}

	res <- ifelse(res <= 0, NA, res)
	
	return(res)
}
