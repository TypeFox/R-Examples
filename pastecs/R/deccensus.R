"deccensus" <-
function(x, type="multiplicative", trend=FALSE) {			# Only the multiplicative model is allowed. Use loess for an additive seasonal decomposition
	# But here we also offer the possibility of using an additive model
	call <- match.call()
	x <- as.ts(x)
	if (is.matrix(x) && ncol(x) != 1)
		stop("only univariate series are allowed")
	# Check the type argument
	TYPES <- c("additive", "multiplicative")
		typeindex <- pmatch(type, TYPES)
		if (is.na(typeindex)) 
			stop("invalid type value")
		if (typeindex == -1) 
			stop("ambiguous type value")
		# make sure type is fully spelled
		type <- switch(typeindex,
				"additive"="additive",
				"multiplicative"="multiplicative")
	if (type == "additive")
		stop("'census' method allows only a multiplicative model. Use 'loess' instead.")
	# create our own specs component
	specs <- list(method="census", type=type)
	# we recuperate units from x
	units <- attr(x, "units")
	# perform filtering
	n <- length(x)
	period <- frequency(x)
    if (period < 2 || n < 3 * period) 
        stop("series is not periodic or has less than three periods")
   # if (period < 2 || n < 4 * period) 
   #     stop("series is not periodic or has less than four periods")
    if (n %% period != 0)
    	stop("the series must have an integer number of periods")
    if(sum(is.na(x)) > 0)
		stop("no missing value allowed for this method")
	
	X <- as.vector(x)
	# 1) This is to test if there is a seasonal fluctuation
	test <- c(NA, 200*X[2:(n-1)]/(X[1:(n-2)]+X[3:n]), NA)
	test <- split(test, cycle(x))
	test.s <- NULL
	for (i in 1:period)
		test.s[i] <- mean(test[[i]], na.rm=TRUE)
	# If any of the values in test.s is NaN or Inf, impossible to apply the CENSUS II method (sparse data?)
	if (sum(is.na(test.s)) > 0 | sum(test.s == Inf) > 0) {
		warning("The 'census' decomposition is not calculable for one series")
		series <- x
		#dimnames(series)[[2]] <- "serie"		# No decomposition was performed
		Amp <- NULL
	} else {
		# 2) A mobile mean with a window equal to frequency is applied
		if (period %% 2 == 0) { 
			filter1 <- c(1/(2*period), rep(1/period, period-1), 1/(2*period)) 
		} else {
			filter1 <- rep(1/period, period)
		}
		A <- filter(x, filter1, method="convolution")

		# 3) First estimation of S x I (except for 6 terms at the beginning and at the end)
		B <- x / A 			#* 100
		B2 <- as.vector(B)
		nc <- n/period
		dim(B2) <- c(period, nc)
	
		# 4) Regulation of B (order 2 mobile mean accross years), and we add two values at each size first
		S2 <- B2
		nr2 <- period/2
		Sa <- S2[1:nr2, 2:nc]
		Sa.lead <- apply(Sa[,1:2], 1, mean)
		Sa.trail <- apply(Sa[,(nc-2):(nc-1)], 1, mean)
		Sa <- cbind(Sa.lead, Sa.lead, Sa, Sa.trail, Sa.trail)
		Sa <- t(apply(Sa, 1, filter, rep(1/5, 5), method="convolution"))
		S2[1:nr2, 2:nc] <- Sa[,3:(nc+1)]
		Sb <- S2[(nr2+1):period, 1:(nc-1)]
		Sb.lead <- apply(Sb[,1:2], 1, mean)
		Sb.trail <- apply(Sb[,(nc-2):(nc-1)], 1, mean)
		Sb <- cbind(Sb.lead, Sb.lead, Sb, Sb.trail, Sb.trail)
		Sb <- t(apply(Sb, 1, filter, rep(1/5, 5), method="convolution"))
		S2[(nr2+1):period, 1:(nc-1)] <- Sb[,3:(nc+1)]
		S <- S2
		dim(S) <- NULL
		# Detection of extreme values:
		s <- sqrt(apply((B2 - S2)^2, 1, sum, na.rm=TRUE)/(nc-2))
		s2 <- rep(1.96*s, nc)
		dim(s2) <- c(period, nc)
		is.xtreme <- abs(B2 - S2) > s2
		is.xtreme[is.na(is.xtreme)] <- FALSE
		# Replacement of extreme values
		B2.rep <- B2
		B2.rep[, 2:(nc-1)] <- (B2[, 1:(nc-2)] + B2[, 2:(nc-1)] + B2[, 3:nc]) / 3
		B2.rep[1:nr2, 2] <- B2.rep[1:nr2, 3]
		B2.rep[(nr2+1):period, 1] <- B2.rep[(nr2+1):period, 2]
		B2.rep[1:nr2, nc] <- B2.rep[1:nr2, nc-1]
		B2.rep[(nr2+1):period, nc-1] <- B2.rep[(nr2+1):period, nc-2]
		B2[is.xtreme] <- B2.rep[is.xtreme]
		# Replace missing values by values for same month and previous/next year
		B2[1:nr2, 1] <- B2[1:nr2, 2]
		B2[(nr2+1):period, nc] <- B2[(nr2+1):period, nc-1]
		# Centering of values
		B2.cent <- rep(apply(B2, 2, sum), period)
		dim(B2.cent) <- c(nc, period)
		B2.cent <- t(B2.cent)
		S1 <- B2 * (period * 100) / B2.cent
				
		# 5) First approximation of S1 with a mobile mean
		S1.lead <- apply(S1[,1:2], 1, mean)
		S1.trail <- apply(S1[,(nc-1):nc], 1, mean)
		S1 <- cbind(S1.lead, S1.lead, S1, S1.trail, S1.trail)
		S1 <- t(apply(S1, 1, filter, c(1/9, 2/9, 3/9, 2/9, 1/9), method="convolution"))
		S1 <- S1[, 3:(nc+2)]
	
		# 6) First approximation of CI1
		dim(S1) <- NULL
		CI1 <- x / S1
				
		# 7) First approximation of C1 (add 7 values at begin and end and apply Spencer filter)
		C1 <- c(rep(mean(CI1[1:4]), 7), CI1, rep(mean(CI1[(n-3):n]), 7))
		C1 <- filter(C1, c(-3,-6,-5,3,21,46,67,74,67,46,21,3,-5,-6,-3)/320, method="convolution")
		C1 <- as.vector(C1[8:(n+7)])
	
		# 8) First approximation of I1
		I1 <- CI1 / C1*100		#*100
		Amp <- sum(abs(I1[2:n] - I1[1:(n-1)]))/(n-1)
		if (is.na(Amp)) {
				warning("The 'census' decomposition is not calculable for one series (Amp == NA)")
				series <- x
		} else {
		
			# 9) Second approximation of SI
			D <- x / C1*100			#*100
			B2 <- as.vector(D)
			dim(B2) <- c(period, nc)
			S2.lead <- apply(B2[,1:2], 1, mean)
			S2.trail <- apply(B2[,(nc-1):nc], 1, mean)
			S2 <- cbind(S2.lead, S2.lead, B2, S2.trail, S2.trail)
			S2 <- t(apply(S2, 1, filter, rep(1/5, 5), method="convolution"))
			S2 <- S2[, 3:(nc+2)]
			# Detection of extreme values:
			s <- sqrt(apply((B2 - S2)^2, 1, sum)/(nc-2))
			s2 <- rep(1.96*s, nc)
			dim(s2) <- c(period, nc)
			is.xtreme <- abs(B2 - S2) > s2
			is.xtreme[is.na(is.xtreme)] <- FALSE
			# Replacement of extreme values
			B2.rep <- B2
			B2.rep[, 2:(nc-1)] <- (B2[, 1:(nc-2)] + B2[, 2:(nc-1)] + B2[, 3:nc]) / 3
			B2.rep[, 1] <- B2.rep[, 2]
			B2.rep[, nc] <- B2.rep[, nc-1]
			B2[is.xtreme] <- B2.rep[is.xtreme]
			# Centering of values
			B2.cent <- rep(apply(B2, 2, sum), period)
			dim(B2.cent) <- c(nc, period)
			B2.cent <- t(B2.cent)
			S1 <- B2 * (period * 100) / B2.cent
		
			# 10) Final season indices
			if (Amp > 2) {
				S1.lead <- apply(S1[,1:3], 1, mean)
				S1.trail <- apply(S1[,(nc-2):nc], 1, mean)
				S1 <- cbind(S1.lead, S1.lead, S1.lead, S1, S1.trail, S1.trail, S1.trail)
				S1 <- t(apply(S1, 1, filter, c(1/16, 2/16, 3/16, 4/16, 3/16, 2/16, 1/16), method="convolution"))
				S <- S1[, 4:(nc+3)]
			} else {
				S1.lead <- apply(S1[,1:2], 1, mean)
				S1.trail <- apply(S1[,(nc-1):nc], 1, mean)
				S1 <- cbind(S1.lead, S1.lead, S1, S1.trail, S1.trail)
				S1 <- t(apply(S1, 1, filter, c(1/9, 2/9, 3/9, 2/9, 1/9), method="convolution"))
				S <- S1[, 3:(nc+2)]
			}
			
			# 11) Final deseasoned series CI
			dim(S) <- NULL
			CI <- x / S	*100
		
			# 12) Final cyclic trend component C
			C <- c(rep(mean(CI[1:4]), 7), CI, rep(mean(CI[(n-3):n]), 7))
			C <- filter(C, c(-3,-6,-5,3,21,46,67,74,67,46,21,3,-5,-6,-3)/320, method="convolution")
			C <- as.vector(C[8:(n+7)])
	
			# 13) Final random component I
			I <- CI / C	*100
			Amp <- sum(abs(I1[2:n] - I1[1:(n-1)]))/(n-1)

			# Concatenate series
			if (trend == FALSE) {
				S <- as.ts(S)
				tsp(S) <- tsp(CI)
				series <- ts.union(CI, S/100)
				dimnames(series)[[2]] <- c("deseasoned", "seasonal")
			} else {
				S <- as.ts(S)
				tsp(S) <- tsp(I)
				C <- as.ts(C)
				tsp(C) <- tsp(I)
				series <- ts.union(C, S/100, I/100)
				dimnames(series)[[2]] <- c("trend", "seasonal", "residuals")
			}
		}
	}
	# create our own 'tsd' structure
	res <- list(ts="series", series=series, test.seasons=test.s, amplitude=Amp, units=units, model.type="multiplicative", specs=specs, call=call)
	class(res) <- "tsd"		# change the class of the object to 'tsd'
	res
}
