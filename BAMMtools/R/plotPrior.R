plotPrior <- function(mcmc, expectedNumberOfShifts = 1, burnin = 0.15, priorCol = 'light blue', postCol = 'red', legendPos = 'topright', ...) {
	
	if (!class(mcmc) %in% c('character', 'data.frame', 'matrix')) {
		stop('mcmc must be either a dataframe or the path to the mcmc_out file.')
	}
	
	if (is.character(mcmc)) {
		mcmc <- read.csv(mcmc, stringsAsFactors = FALSE)
	}
	
	#drop burnin
	mcmc2 <- mcmc[floor(burnin * nrow(mcmc)):nrow(mcmc),]
	
	#get prior distribution of shifts
	obsK <- seq(from = 0, to = max(mcmc2[,"N_shifts"]), by = 1)
	prior <- sapply(obsK, prob.k, poissonRatePrior = 1/expectedNumberOfShifts)
	prior <- data.frame(N_shifts = obsK, prob = prior)
	
	#get posterior distribution of shifts
	posterior <- sapply(obsK, function(x) length(which(mcmc2[,'N_shifts'] == x))) / nrow(mcmc2)
	names(posterior) <- obsK
	posterior <- data.frame(N_shifts = names(posterior), prob = posterior)

	barplot(prior[,2], names.arg = prior[,1], ylim = c(0, max(c(prior[,2], posterior[,2]))), border = 'black', col = priorCol, xlab = 'n shifts', ...)
	barplot(posterior[,2], add = TRUE, border = 'black', col = BAMMtools::transparentColor(postCol, 0.4), axes=FALSE)
	legend(x = legendPos, y = NULL, legend = c('prior','posterior'), fill = c(priorCol, BAMMtools::transparentColor(postCol, 0.4)), bty = 'n', cex=1.5)
	
	invisible(cbind(N_shifts = prior$N_shifts, priorProbs = prior$prob, postProbs = posterior$prob))	
}


prob.k <- function(k, poissonRatePrior) {
	Denom <- (poissonRatePrior + 1) ^ (k + 1)
	Prob <- poissonRatePrior / Denom
	return(Prob)
}
