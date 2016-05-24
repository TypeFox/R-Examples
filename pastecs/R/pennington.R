"pennington" <-
function(x, calc="all", na.rm=FALSE) {
	# Calculation of Pennington's mean
	pen.mean <- function(n, m, MeanLogVal, VarLogVal) {
		# Calculation of Gm(t)
		Gmt <- function(m, t) {
			# This function calculate a single j term for the sum
			SumTerm <- function(m, t, j) {ST <- exp((2*j-1)*log(m-1)+j*log(t)-j*log(m)-sum(log(c(2:j-1)*2-1+m))-sum(log(c(1:j)))); ST}
			# Sum is calculated iteratively to a precision of 0.01%
			# If t=0, Gmt=1. m must not be 0 or 1, but this is treated in the body of the main function
			if (t==0) Sum <- 1 else {
				Sum <- 1+t*(m-1)/m	# we already put the first terms in Sum
				it <- 1		# iteration count (will start at 2)
				del <- 1	# iterative adjustment
				while (del>0.0001 && (it <- it+1)<20) {
					Sum2 <- Sum+SumTerm(m, t, it)
					del <- abs((Sum-Sum2)/Sum)
					Sum <- Sum2
				}
			}
			Sum
		}
		
		# No non-zero values => return 0
		if (m==0) PenMean <- 0 else {
			# Only one non-zero value => return x1/n
			if (m==1) PenMean <- exp(MeanLogVal)/n else {
				# Calculation with the full formula
				PenMean <- m/n*exp(MeanLogVal)*Gmt(m,VarLogVal/2)
			}
		}
		PenMean
	}
	
	# Calculation of Pennington's variance
	pen.var <- function(n, m, MeanLogVal, VarLogVal) {
		# Calculation of Gm(t)
		Gmt <- function(m, t) {
			# This function calculate a single j term for the sum
			SumTerm <- function(m, t, j) {ST <- exp((2*j-1)*log(m-1)+j*log(t)-j*log(m)-sum(log(c(2:j-1)*2-1+m))-sum(log(c(1:j)))); ST}
			# Sum is calculated iteratively to a precision of 0.01%
			# If t=0, Gmt=1. m must not be 0 or 1, but this is treated in the body of the main function
			if (t==0) Sum <- 1 else {
				Sum <- 1+t*(m-1)/m	# we already put the first terms in Sum
				it <- 1		# iteration count (will start at 2)
				del <- 1	# iterative adjustment
				while (del>0.0001 && (it <- it+1)<20) {
					Sum2 <- Sum+SumTerm(m, t, it)
					del <- abs((Sum-Sum2)/Sum)
					Sum <- Sum2
					# cat(it, Sum, "\n")	# To get all iterations
				}
			}
			Sum
		}

		# No non-zero values => return 0
		if (m==0) PenVar <- 0 else {
			# Only one non-zero value => return x1^2/n
			if (m==1) PenVar <- (exp(MeanLogVal))^2/n else {
				# Calculation with the full formula
				PenVar <- m/n*exp(2*MeanLogVal)*(Gmt(m,2*VarLogVal)-(m-1)/(n-1)*Gmt(m,(m-2)/(m-1)*VarLogVal))
			}
		}
		PenVar
	}
		
	# Calculation of Pennington's variance of the mean
	pen.mean.var <- function(n, m, MeanLogVal, VarLogVal) {
		# Calculation of Gm(t)
		Gmt <- function(m, t) {
			# This function calculate a single j term for the sum
			SumTerm <- function(m, t, j) {ST <- exp((2*j-1)*log(m-1)+j*log(t)-j*log(m)-sum(log(c(2:j-1)*2-1+m))-sum(log(c(1:j)))); ST}
			# Sum is calculated iteratively to a precision of 0.01%
			# If t=0, Gmt=1. m must not be 0 or 1, but this is treated in the body of the main function
			if (t==0) Sum <- 1 else {
				Sum <- 1+t*(m-1)/m	# we already put the first terms in Sum
				it <- 1		# iteration count (will start at 2)
				del <- 1	# iterative adjustment
				while (del>0.0001 && (it <- it+1)<20) {
					Sum2 <- Sum+SumTerm(m, t, it)
					del <- abs((Sum-Sum2)/Sum)
					Sum <- Sum2
					# cat(it, Sum, "\n")	# To get all iterations
				}
			}
			Sum
		}

		# No non-zero values => return 0
		if (m==0) PenMeanVar <- 0 else {
			# Only one non-zero value => return x1^2/n
			if (m==1) PenMeanVar <- (exp(MeanLogVal))^2/n else {
				# Calculation with the full formula
				PenMeanVar <- m/n*exp(2*MeanLogVal)*(m/n*(Gmt(m,VarLogVal/2))^2-(m-1)/(n-1)*Gmt(m,(m-2)/(m-1)*VarLogVal))
			}
		}
		PenMeanVar
	}
		
	
	# This is the core part of pennington!
	# If na.rm=FALSE and some missing data, must return NA
	if (na.rm==FALSE & length(x[!is.na(x)])<length(x)) IsNa <- TRUE else IsNa <- FALSE
	# N is the nbr of values (excluding missing ones)
	N <- sum(as.numeric(!is.na(x)))
	# If no remaining values after eliminating missing values, return NA
	if (N==0) IsNa <- TRUE
	# End of tests => either return NA or further calculate
	if (IsNa==TRUE) {
		if (calc=="all") Result <- unlist(list(mean=NA, var=NA, mean.var=NA)) else Result <- NA
	} else {		# Calculate
		# LogVal contains log of all non-zero and non-missing values
		LogVal <- log(x[!is.na(x) & x>0])
		# M is the number of non-zero values
		M <- length(LogVal)
		# Calculation of mean and variance
		MeanLogVal <- mean(LogVal)
		if (M == 0) VarLogVal <- NaN else VarLogVal <- var(LogVal)
		# Depending on the calc arg, we calculate the mean, var or var.mean
		Result <- switch(calc,
			mean = pen.mean(N, M, MeanLogVal, VarLogVal),
			var = pen.var(N, M, MeanLogVal, VarLogVal),
			mean.var = pen.mean.var(N, M, MeanLogVal, VarLogVal),
			all = unlist(list(mean=pen.mean(N, M, MeanLogVal, VarLogVal), var=pen.var(N, M, MeanLogVal, VarLogVal), mean.var=pen.mean.var(N, M, MeanLogVal, VarLogVal))))
	}
	Result	
}
