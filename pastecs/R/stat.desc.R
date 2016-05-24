"stat.desc" <-
function(x, basic=TRUE, desc=TRUE, norm=FALSE, p=.95) {
	# This function performs all calculations on a single vector
	# Missing data allowed and stripped out before calculations
	stat.desc.vec <- function(x, basic, desc, norm, p) {
		# If x is a	list, we transform it into a vector
		x <- unlist(x)
		if (!is.numeric(x)) {	# Not a numeric vector!
			Nbrval <- NA; Nbrnull <- NA; Nbrna <-NA
			Median <- NA; Mean <- NA; StdDev <- NA
			if (basic==TRUE) {
				Res1 <- list(nbr.val=NA, nbr.null=NA, nbr.na=NA, min=NA, max=NA, range=NA, sum=NA)
			} else Res1 <- NULL
			if (desc==TRUE) {
				CIMean <- NA; names(CIMean) <- p		# To get CI.mean.pValue
				Res2 <- list(median=NA, mean=NA, SE.mean=NA, CI.mean=NA, var=NA, std.dev=NA, coef.var=NA)
			} else Res2 <- NULL
			if (norm==TRUE) {
				Res3 <- list(skewness=NA, skew.2SE=NA, kurtosis=NA, kurt.2SE=NA, normtest.W=NA, normtest.p=NA)
			} else Res3 <- NULL
		} else {			# Vector contains numbers, we can perform calcs
			Nbrna <- sum(as.numeric(is.na(x)))
			# We could use na.rm=TRUE everywhere, but it is faster
			# to remove all missing values once at the beginning
			x <- x[!is.na(x)]
			Nbrval <- length(x)
			Nbrnull <- sum(as.numeric(x==0))
			if (basic==TRUE) {
				Min <- min(x)
				Max <- max(x)
				Range <- Max-Min
				Sum <- sum(x)
				Res1 <- list(nbr.val=Nbrval, nbr.null=Nbrnull, nbr.na=Nbrna, min=Min, max=Max, range=Range, sum=Sum)
			} else Res1 <- NULL
			Median <- median(x); names(Median) <- NULL	# To correct a bug!?
			Mean <- mean(x)
			Var <- var(x)
			StdDev <- sqrt(Var)
			SEMean <- StdDev/sqrt(Nbrval)
			if (desc==TRUE) {
				CIMean <- qt((0.5+p/2), (Nbrval-1))*SEMean
				# (0.5+p/2) because we need a two-sided t distribution for the CI!!!
				names(CIMean) <- p
				CoefVar <- StdDev/Mean
				Res2 <- list(median=Median, mean=Mean, SE.mean=SEMean, CI.mean=CIMean, var=Var, std.dev=StdDev, coef.var=CoefVar)
			} else Res2 <- NULL
			if (norm==TRUE) {
				Skew <- sum((x-mean(x))^3)/(length(x)*sqrt(var(x))^3)			# From e1071
				Kurt <- sum((x-mean(x))^4)/(length(x)*var(x)^2) - 3				# Idem
				SE <- sqrt(6*Nbrval*(Nbrval-1)/(Nbrval-2)/(Nbrval+1)/(Nbrval+3))
				# +/- sqrt(6/Nbrval) for Nbrval>100, see Sokal & Rohlf, p. 139
				Skew.2SE <- Skew/(2*SE)
				# If abs(Skew.2SE)>1 then Skew is signific., see Systat 9 p. I-212
				SE <- sqrt(24*Nbrval*((Nbrval-1)^2)/(Nbrval-3)/(Nbrval-2)/(Nbrval+3)/(Nbrval+5))
				# +/- sqrt(24/Nbrval) for Nbrval > 150, same ref.
				Kurt.2SE <- Kurt/(2*SE)
				# Same remark as for Skew.2SE!
				# This is the Shapiro-Wilk test of normality
				Ntest <- shapiro.test(x)
				Ntest.W <- Ntest$statistic; names(Ntest.W) <- NULL
				Ntest.p <- Ntest$p.value
				Res3 <- list(skewness=Skew, skew.2SE=Skew.2SE, kurtosis=Kurt, kurt.2SE=Kurt.2SE, normtest.W=Ntest.W, normtest.p=Ntest.p)
			} else Res3 <- NULL
		}
		# We collect all results together
		Res <- unlist(c(Res1, Res2, Res3))
		# If basic, desc & norm are all false, we return just minimal calculations
		if (length(Res)==0) Res <- unlist(list(nbr.val=Nbrval, nbr.null=Nbrnull, nbr.na=Nbrna, median=Median, mean=Mean, std.dev=StdDev))
		Res
	}
	
	# This is the body of stat.desc
	Basic <- basic; Desc <- desc; Norm <- norm; P <- p
	# If x is a vector, stat.desc returns a vector with results
	# TO DO: allow also time series (ts) => is.data.frame, or stat.desc.vec?
	if (is.vector(x)) stat.desc.vec(x, Basic, Desc, Norm, P) else {
		# If x is not a vector, it is treated as a data frame
		# A result will be returned in a data frame with corresponding columns
		x <- as.data.frame(x)
		# We keep the same column headers
		NamesV <- names(x)
		StatM <- NULL
		# Calculation is performed alternatively on each column
		for (i in 1:ncol(x)) {
			StatV <- stat.desc.vec(x[i], Basic, Desc, Norm, P)
			# The next if condition to avoid error at the first step
			if (is.null(StatM)==TRUE) StatM <- data.frame(StatV) else StatM <- cbind(StatM, StatV)
		}
		# We change names of columns to match the original data frame
		names(StatM) <- NamesV
		StatM
	}
}
