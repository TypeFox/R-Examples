
emInit <- function(y, g, conds, lib.size, lib.type, alg.type = "EM", 
	init.runs, init.iter, fixed.lambda, equal.proportions, 
	verbose, s=NA) {

	if(is.matrix(y) == FALSE & is.data.frame(y) == FALSE) 
		stop(paste(sQuote("y"), "must be a matrix"))
	if(min(y) < 0)
		stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
	if(sum(round(y)) != sum(y)) 
		stop(paste(sQuote("y"), "must be a matrix made up of nonnegative counts"))
	if(min(rowSums(y)) == 0)
		stop(paste("at least one observation in", sQuote("y"), "contains all 0's and must be removed from the data"))
	if(length(g) != 1)
		stop(paste(sQuote("g"), "(the number of clusters) must be a nonnegative integer"))
	if(g <= 0 | round(g) != g) 
		stop(paste(sQuote("g"), "(the number of clusters) must be a nonnegative integer"))
	if(is.vector(conds) == FALSE | length(conds) != ncol(y))
		stop(paste(sQuote("conds"), "must be a vector the same length as the number of columns in", sQuote ("y")))
	if(is.logical(lib.size) == FALSE)
		stop(paste(sQuote("libsize"), "must be", dQuote("TRUE"), "(PMM-II) or", 
			dQuote("FALSE"), "(PMM-I)"))
	if(alg.type != "EM" & alg.type != "CEM")
		stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
	if(length(alg.type) > 1)
		stop(paste(sQuote("alg.type"), "must be one of", dQuote("EM"), "or", dQuote("CEM")))
	if(init.runs < 1 | length(init.runs) > 1 | round(init.runs) != init.runs) 
		stop(paste(sQuote("init.runs"), "must be a positive integer"))
	if(is.logical(verbose) == FALSE)
		stop(paste(sQuote("verbose"), "must be", dQuote("TRUE"), "or", dQuote("FALSE")))
	init.type1 <- "kmeans"
	lambda.init.all <- vector("list", init.runs)
	pi.init.all <- vector("list", init.runs)
	criterion.all <- rep(NA, init.runs)
	for(start in 1:init.runs) {
		em.init <- PoisMixClus(y = y, g = g, lib.size = lib.size, 
			lib.type = lib.type, conds = conds, 
			init.type = init.type1, alg.type = alg.type, iter = init.iter,
			fixed.lambda = fixed.lambda, equal.proportions = equal.proportions, s = s,
			wrapper=FALSE)
			lambda.init.all[[start]] <- em.init$lambda
			pi.init.all[[start]] <- em.init$pi
			criterion.all[start] <- em.init$log.like
		if(verbose == TRUE) print(paste("Initialization:", start))
	}
	## If all criterion values are equal to NaN, then arbitrarily choose the first one
	if(sum(is.na(criterion.all)) == length(criterion.all)) {
		final.choice <- 1;
	}
	## If two of the criterion values are exactly the same, pick only the first
	if(sum(is.na(criterion.all)) != length(criterion.all)) {
		final.choice <- which(criterion.all == min(criterion.all, na.rm = TRUE))[1]
	}
	lambda.init <- lambda.init.all[[final.choice]]
	pi.init <- pi.init.all[[final.choice]]
	return(list(pi.init = pi.init, lambda.init = lambda.init))
}
