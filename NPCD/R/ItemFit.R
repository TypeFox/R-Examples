###############################################################################
# ItemFit:                                                                    #
#                                                                             #
# Define functions to generate item fit statistics for various classes of     #
# outputs generated from the functions in this package, including AlphaNP,    #
# AlphaMLE, and JMLE. Currently supports RMSEA and Chi-square.                #
#                                                                             #
###############################################################################

#install.packages("R.oo")
#library(R.oo)

# inputs: x: the output from other functions
# outputs: item fit statistics and hypothesis testing result

################################################################################

ItemFit <- function(x, model = NULL, par = NULL) {

	# Extract necessary information
	
	Q <- x$Q
	Y <- x$Y

	nitem <- nrow(Q)
	natt <- ncol(Q)
	nperson <- nrow(Y)
	pattern <- AlphaPermute(natt)
	npattern <- nrow(pattern)

	# Different treatments for classes
	
	if (class(x) %in% c("AlphaNP", "AlphaMLE", "JMLE")) 
		est.class <- x$est.class

	if (class(x) %in% c("AlphaMLE", "ParMLE", "JMLE")) 
		model <- x$model

	if (class(x) == "AlphaMLE") 
		par <- x$par

	if (class(x) == "ParMLE") {
		par <- x
		est.class <- rep(NA, nperson)
		for (j in 1:npattern) est.class[apply(x$alpha, 1, function(m) all(m == pattern[j, ]))] <- j
	}

	if (class(x) == "JMLE") 
		par <- x$par.est

	# Different treatments for models
	
	if (model %in% c("DINA", "DINO", "NIDA", "GNIDA")) {
		slip <- par$slip
		guess <- par$guess
	}

	if (model == "RRUM") {
		pi <- par$pi
		r <- par$r
	}

	if (model %in% c("DINA", "DINO")) 
		Chisq.df <- npattern - 2
	if (model %in% c("NIDA", "GNIDA")) 
		Chisq.df <- npattern - natt
	if (model == "RRUM") 
		Chisq.df <- npattern - natt - 1

	# fit indices
	
	class.freq <- rep(NA, npattern)
	for (j in 1:npattern) class.freq[j] <- sum(est.class == j)
	class.prop <- class.freq/sum(class.freq)

	RMSEA <- Chisq.p <- rep(NA, nitem)
	my.Chisq <- rep(0, nitem)

	for (i in 1:nitem) {
		if (model %in% c("DINA", "DINO")) {
			par.tmp <- list(slip = slip[i], guess = guess[i])
		}
		if (model == "GNIDA") {
			par.tmp <- list(slip = slip[i, ], guess = guess[i, ])
		}
		if (model %in% c("NIDA")) {
			par.tmp <- list(slip = slip, guess = guess)
		}
		if (model == "RRUM") {
			par.tmp <- list(pi = pi[i], r = r[i, ])
		}

		RMSEA.tmp <- 0

		for (j in 1:npattern) {

			P <- CDP(Q[i, ], par.tmp, pattern[j, ], model)
			P.obs <- mean(Y[est.class == j, i])

			RMSEA.tmp <- RMSEA.tmp + (P - P.obs)^2 * class.prop[j]

			O <- P.obs * class.freq[j]
			E <- P * class.freq[j]
			N <- class.freq[j]
			my.Chisq[i] <- my.Chisq[i] + N * (O - E)^2/E/(N - E)
		}

		RMSEA[i] <- sqrt(RMSEA.tmp)
		Chisq.p[i] <- 1 - pchisq(my.Chisq[i], Chisq.df)
	}

	out <- data.frame(RMSEA = RMSEA, Chisq = my.Chisq, Chisq.p = Chisq.p, Chisq.df = Chisq.df)

	# Row names and column names of the output
	
	my.rowname <- rep(NA, nitem)
	for (i in 1:nitem) my.rowname[i] <- paste("Item", i)
	rownames(out) <- my.rowname
	colnames(out) <- c("RMSEA", "Chisq", "Chisq p-value", "Chisq df")

	return(out)
}