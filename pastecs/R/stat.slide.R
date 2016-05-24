"stat.slide" <-
function(x, y, xcut=NULL, xmin=min(x), n=NULL, frequency=NULL, deltat=1/frequency, basic=FALSE, desc=FALSE, norm=FALSE, pen=FALSE, p=.95) {
	# This function performs all calculations on a single bloc
	# Missing data allowed and stripped out before calculations
	stat.slide.cut <- function(y, xmin, xmax, basic, desc, norm, pen, p) {
		# If y is a	list, we transform it into a vector
		y <- unlist(y)
		if (!is.numeric(y)) {	# Not a numeric vector!
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
			if (pen==TRUE) {
				Res4 <- list(pos.median=NA, pos.mean=NA, geo.mean=NA, pen.mean=NA, pen.var=NA, pen.std.dev=NA, pen.mean.var=NA)
			} else Res4 <- NULL
		} else {			# Vector contains numbers, we can perform calcs
			# First two entries are xmin and xmax for each cut
			Res0 <- list(xmin=xmin, xmax=xmax)
			Nbrna <- sum(as.numeric(is.na(y)))
			# We could use na.rm=TRUE everywhere, but it is faster
			# to remove all missing values once at the beginning
			y <- y[!is.na(y)]
			Nbrval <- length(y)
			Nbrnull <- sum(as.numeric(y==0))
			if (basic==TRUE) {
				Min <- min(y)
				Max <- max(y)
				Range <- Max-Min
				Sum <- sum(y)
				Res1 <- list(nbr.val=Nbrval, nbr.null=Nbrnull, nbr.na=Nbrna, min=Min, max=Max, range=Range, sum=Sum)
			} else Res1 <- NULL
			Median <- median(y); names(Median) <- NULL	# To correct a bug!?
			Mean <- mean(y)
			Var <- var(y)
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
				Skew <- sum((y-mean(y))^3)/(length(y)*sqrt(var(y))^3)			# From e1071
				Kurt <- sum((y-mean(y))^4)/(length(y)*var(y)^2) - 3				# Idem
				SE <- sqrt(6*Nbrval*(Nbrval-1)/(Nbrval-2)/(Nbrval+1)/(Nbrval+3))
				# +/- sqrt(6/Nbrval) for Nbrval>100, see Sokal & Rohlf, p. 139
				Skew.2SE <- Skew/(2*SE)
				# If abs(Skew.2SE)>1 then Skew is signific., see Systat 9 p. I-212
				SE <- sqrt(24*Nbrval*((Nbrval-1)^2)/(Nbrval-3)/(Nbrval-2)/(Nbrval+3)/(Nbrval+5))
				# +/- sqrt(24/Nbrval) for Nbrval > 150, same ref.
				Kurt.2SE <- Kurt/(2*SE)
				# Same remark as for Skew.2SE!
				# This is the Shapiro-Wilk test of normality
				# Now done with Depends: field require(stats)		# For Shapiro-Wilk normality tests in R
				Ntest <- shapiro.test(y)
				Ntest.W <- Ntest$statistic; names(Ntest.W) <- NULL
				Ntest.p <- Ntest$p.value
				Res3 <- list(skewness=Skew, skew.2SE=Skew.2SE, kurtosis=Kurt, kurt.2SE=Kurt.2SE, normtest.W=Ntest.W, normtest.p=Ntest.p)
			} else Res3 <- NULL
			if (pen==TRUE) {
				ypos <- y[y>0]
				PosMedian <- median(ypos); names(PosMedian) <- NULL
				PosMean <- mean(ypos)
				GeoMean <- exp(mean(log(ypos)))
				Pen <- pennington(y, calc="all")
				names(Pen) <- NULL
				PMean <- Pen[1]
				PVar <- Pen[2]
				PStdDev <- sqrt(PVar)
				PMeanVar <- Pen[3]
				Res4 <- list(pos.median=PosMedian, pos.mean=PosMean, geo.mean=GeoMean, pen.mean=PMean, pen.var=PVar, pen.std.dev=PStdDev, pen.mean.var=PMeanVar)
			} else Res4 <- NULL
		}
		# We collect all results together
		Res <- unlist(c(Res0, Res1, Res2, Res3, Res4))
		# If basic, desc, norm & pen are all false, we return just minimal calculations
		if (length(Res)==2) Res <- unlist(list(xmin=xmin, xmax=xmax, nbr.val=Nbrval, nbr.null=Nbrnull, nbr.na=Nbrna, min=min(y), max=max(y), median=Median, mean=Mean, std.dev=StdDev))
		Res
	}
	
	# This is the body of stat.slide
	call <- match.call()
	Basic <- basic; Desc <- desc; Norm <- norm; Pen <- pen; P <- p
	x <- xy.coords(x, y)
	y <- x$y
	x <- x$x
	# Make sure data are sorted in increasing order according to x
	srt <- sort.list(x)
	x <- x[srt]
	y <- y[srt]
	# Check arguments
	if (!is.numeric(x) || !is.numeric(y)) 
		stop("stat.slide: x and y must be numeric")
	nx <- length(x)
	if (nx != length(y)) 
	    stop("x and y must have equal lengths")
	if (nx < 3) 
	    stop("stat.slide requires at least three values for x and y")
	# Eliminate entries with missing values
	ok <- !(is.na(x) | is.na(y))
	x <- x[ok]
	y <- y[ok]
	nx <- length(x)
	if (nx < 3) 
	    stop("stat.slide requires at least three non-missing values for x and y")
	# The different sequences are provided in xcut, or are calculated from xmin, frequency/deltat and n
	if (is.null(xcut)) {									# calculate regular sequences
		if (is.null(deltat) | deltat <= 0) 
	    	stop("stat.slide requires deltat > 0")
		# if n is not defined, then it is calculated
		if (is.null(n)) {
			n <- ceiling((max(x)-min(x))/deltat)
		}
		xcut <- 0:n * deltat + xmin							# vector of xcut regular sequence
	}
	StatM <- NULL
	ncut <- length(xcut) - 1
	# Calculation is performed alternatively on each bloc
	for (i in 1:ncut) {
		StatV <- stat.slide.cut(y[x >= xcut[i] & x < xcut[i+1]], xcut[i], xcut[i+1], Basic, Desc, Norm, Pen, P)
		# The next if condition to avoid error at the first step
		if (is.null(StatM)==TRUE) StatM <- data.frame(StatV) else StatM <- cbind(StatM, StatV)
	}
	# We change names of columns to represent the cuts
	NamesV <- paste("[", xcut[1:ncut], ",", xcut[2:(ncut+1)], "[", sep="")
	names(StatM) <- NamesV
	StatM
	res <- list(stat=StatM, x=x, y=y, xcut=xcut, call=call)			# Create a list containing the result
	class(res) <- "stat.slide"								# and turn it into a 'stat.slide' object
	res
}
