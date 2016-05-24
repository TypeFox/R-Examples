## /*****************************************************************************
##	* SIENA: Simulation Investigation for Empirical Network Analysis
##	*
##	* Web: http://www.stats.ox.ac.uk/~snijders/siena
##	*
##	* File: sienaGOF.r
##	*
##	* Description: This file contains the code to assess goodness of fit:
##	* the main function sienaGOF, the plot method,
##	* and auxiliary statistics and extractor functions.
##	* Written by Josh Lospinoso, modifications by Tom Snijders.
##	*
##	****************************************************************************/

##@sienaGOF siena07 Does test for goodness of fit
sienaGOF <- function(
		sienaFitObject,	auxiliaryFunction,
		period=NULL, verbose=FALSE, join=TRUE, twoTailed=FALSE,
		cluster=NULL, robust=FALSE,
		groupName="Data1", varName, ...)
	{
	require(MASS)
	## require(Matrix)
	##	Check input
	if (! sienaFitObject$returnDeps)
	{
		stop("You must instruct siena07 to return the simulated networks")
	}
	iterations <- length(sienaFitObject$sims)
	if (iterations < 1)
	{
		stop("You need at least one iteration.")
	}
	if (missing(varName))
	{
		stop("You need to supply the parameter <<varName>>.")
	}
	if (missing(auxiliaryFunction))
	{
		stop("You need to supply the parameter <<auxiliaryFunction>>.")
	}
	groups <- length(sienaFitObject$f$groupNames)
	if (verbose)
	{
		if (groups <= 1)
		{
			cat("Detected", iterations, "iterations and", groups, "group.\n")
		}
		else
		{
			cat("Detected", iterations, "iterations and", groups, "groups.\n")
		}
	}

	if (is.null(period) )
	{
		period <- 1:(attr(sienaFitObject$f[[1]]$depvars[[1]], "netdims")[3] - 1)
	}

	 obsStatsByPeriod <- lapply(period, function (j) {
						matrix(
						auxiliaryFunction(NULL,
								sienaFitObject$f,
				sienaFitObject$sims, j, groupName, varName, ...)
						, nrow=1)
				})
	if (join)
	{
		obsStats <- Reduce("+", obsStatsByPeriod)
		obsStats <- list(Joint=obsStats)
	}
	else
	{
		obsStats <- obsStatsByPeriod
		names(obsStats) <- paste("Period", period)
	}
	plotKey <- names(auxiliaryFunction(NULL, sienaFitObject$f,
				sienaFitObject$sims, 1, groupName, varName, ...))
	class(obsStats) <- "observedAuxiliaryStatistics"
	attr(obsStats,"auxiliaryStatisticName") <-
			deparse(substitute(auxiliaryFunction))
	attr(obsStats,"joint") <- join

	##	Calculate the simulated auxiliary statistics
	if (verbose)
	{
		if (length(period) <= 1)
		{
			cat("Calculating auxiliary statistics for period", period, ".\n")
		}
		else
		{
			cat("Calculating auxiliary statistics for periods", period, ".\n")
		}
	}

	if (!is.null(cluster)) {
		ttcSimulation <- system.time(simStatsByPeriod <- lapply(period,
			function (j) {
				simStatsByPeriod <- parSapply(cluster, 1:iterations,
					function (i)
						{ auxiliaryFunction(i,
									sienaFitObject$f,
									sienaFitObject$sims, j, groupName, varName, ...)
							if (verbose && (i %% 100 == 0) )
								{
								cat("  > Completed ", i," calculations\n")
								flush.console()
								}
						})
							simStatsByPeriod <- matrix(simStatsByPeriod,
								ncol=iterations)
							dimnames(simStatsByPeriod)[[2]] <-	1:iterations
							t(simStatsByPeriod)
						}))
	}
	else
	{
		ttcSimulation <- system.time( simStatsByPeriod <- lapply(period,
					function (j) {
						simStatsByPeriod <- sapply(1:iterations, function (i)
						{
							if (verbose && (i %% 100 == 0) )
								{
								cat("  > Completed ", i,
										" calculations\n")
								flush.console()
								}
								auxiliaryFunction(i,
										sienaFitObject$f,
										sienaFitObject$sims, j, groupName, varName, ...)
						})
					simStatsByPeriod <-
							matrix(simStatsByPeriod, ncol=iterations)
					dimnames(simStatsByPeriod)[[2]] <-	1:iterations
					t(simStatsByPeriod)
					})
	  )
	}

	## Aggregate by period if necessary to produce simStats
	if (join)
	{
		simStats <- Reduce("+", simStatsByPeriod)
		simStats <- list(Joint=simStats)
	}
	else
	{
		simStats <- simStatsByPeriod
		names(simStats) <- paste("Period",period)
	}
	class(simStats) <- "simulatedAuxiliaryStatistics"
	attr(simStats,"auxiliaryStatisticName") <-
			deparse(substitute(auxiliaryFunction))
	attr(simStats,"joint") <- join
	attr(simStats,"time") <- ttcSimulation

	applyTest <-  function (observed, simulated)
	{
		if (class(simulated) != "matrix")
		{
			stop("Invalid input.")
		}
		if (class(observed) != "matrix")
		{
			observed <- matrix(observed,nrow=1)
		}
		if (class(observed) != "matrix")
		{
			stop("Observation must be a matrix.")
		}
		if (ncol(observed) != ncol(simulated))
		{
			stop("Dimensionality of function parameters do not match.")
		}
		observations <- nrow(observed)
	#	simulations<-nrow(simulated)
		variates<-ncol(simulated)
		if (robust) {
			a <- cov.rob(simulated)$cov
		}
		else
		{
			a <- cov(simulated)
		}
		ainv <- ginv(a)
		arank <- rankMatrix(a)
		expectation <- colMeans(simulated);
		centeredSimulations <- scale(simulated, scale=FALSE)
		if (variates==1)
		{
			centeredSimulations <- t(centeredSimulations)
		}
		mhd <- function(x)
		{
			x %*% ainv %*% x
		}
		simTestStat <- apply(centeredSimulations, 1, mhd)
		centeredObservations <- observed - expectation
		obsTestStat <- apply(centeredObservations, 1, mhd)
		if (twoTailed)
		{
			p <- sapply(1:observations, function (i)
						1 - abs(1 - 2 * sum(obsTestStat[i] <=
						simTestStat)/length(simTestStat)) )
		}
		else
		{
			p <- sapply(1:observations, function (i)
				sum(obsTestStat[i] <= simTestStat) /length(simTestStat))
		}

		ret <- list( p = p,
				SimulatedTestStat=simTestStat,
				ObservedTestStat=obsTestStat,
				TwoTailed=twoTailed,
				Simulations=simulated,
				Observations=observed,
				InvCovSimStats=a,
				Rank=arank)
		class(ret) <- "sienaGofTest"
		attr(ret,"auxiliaryStatisticName") <-
				attr(obsStats,"auxiliaryStatisticName")
		attr(ret, "key") <- plotKey
		ret
	}

	res <- lapply(1:length(simStats),
					function (i) {
				 applyTest(obsStats[[i]], simStats[[i]]) })
	mhdTemplate <- rep(0, sum(sienaFitObject$test))
	names(mhdTemplate) <- rep(0, sum(sienaFitObject$test))

	JoinedOneStepMHD <- mhdTemplate
	JoinedPartialOneStepMHD <- mhdTemplate
	JoinedGmmMhdValue <- mhdTemplate

	OneStepMHD <- lapply(period, function(i) (mhdTemplate))
	PartialOneStepMHD <- lapply(period, function(i) (mhdTemplate))
	GmmMhdValue <- lapply(period, function(i) (mhdTemplate))

	obsMhd <- NULL

	ExpStat <- lapply(period, function(i) {
				colMeans(simStatsByPeriod[[i]])
			})
	OneStepSpecs <- matrix(0, ncol=sum(sienaFitObject$test),
			nrow=length(sienaFitObject$theta))
	PartialOneStepSpecs <- matrix(0, ncol=sum(sienaFitObject$test),
			nrow=length(sienaFitObject$theta))
	GmmOneStepSpecs <- matrix(0, ncol=sum(sienaFitObject$test),
			nrow=length(sienaFitObject$theta))
	if (robust) {
		covInvByPeriod <- lapply(period, function(i) ginv(
							cov.rob(simStatsByPeriod[[i]]) ))
	}
	else
	{
		covInvByPeriod <- lapply(period, function(i) ginv(
							cov(simStatsByPeriod[[i]]) ))
	}

	obsMhd <- sapply(period, function (i) {
				 (obsStatsByPeriod[[i]] - ExpStat[[i]])	 %*%
						covInvByPeriod[[i]] %*%
						t(obsStatsByPeriod[[i]] - ExpStat[[i]] )
			})

	if (sum(sienaFitObject$test) > 0) {
		effectsObject <- sienaFitObject$requestedEffects
		nSims <- sienaFitObject$Phase3nits
		for (i in period) {
			names(OneStepMHD[[i]]) <-
					effectsObject$effectName[sienaFitObject$test]
			names(PartialOneStepMHD[[i]]) <-
					effectsObject$effectName[sienaFitObject$test]
			names(GmmMhdValue[[i]]) <-
					effectsObject$effectName[sienaFitObject$test]
		}
		names(JoinedOneStepMHD) <-
				effectsObject$effectName[sienaFitObject$test]
		names(JoinedPartialOneStepMHD) <-
				effectsObject$effectName[sienaFitObject$test]
		names(JoinedGmmMhdValue) <-
				effectsObject$effectName[sienaFitObject$test]

		rownames(OneStepSpecs) <- effectsObject$effectName
		colnames(OneStepSpecs) <- effectsObject$effectName[sienaFitObject$test]
		rownames(PartialOneStepSpecs) <- effectsObject$effectName
		colnames(PartialOneStepSpecs) <-
				effectsObject$effectName[sienaFitObject$test]
		rownames(GmmOneStepSpecs) <- effectsObject$effectName
		colnames(GmmOneStepSpecs) <-
				effectsObject$effectName[sienaFitObject$test]

		counterTestEffects <- 0
		for(index in which(sienaFitObject$test)) {
			if (verbose) {
				cat("Estimating test statistic for model including ",
						effectsObject$effectName[index], "\n")
			}
			counterTestEffects <- counterTestEffects + 1
			effectsToInclude <- !sienaFitObject$test
			effectsToInclude[index] <- TRUE
			theta0 <- sienaFitObject$theta
			names(theta0) <- effectsObject$effectName
			theta0 <- theta0[effectsToInclude]
			obsSuffStats <-
					t(sienaFitObject$targets2[effectsToInclude, , drop=FALSE])
			G <- sienaFitObject$sf2[, , effectsToInclude, drop=FALSE] -
					rep(obsSuffStats, each=nSims)
			sigma <- cov(apply(G, c(1, 3), sum))
			SF <- sienaFitObject$ssc[ , , effectsToInclude, drop=FALSE]
			dimnames(SF)[[3]] <- effectsObject$effectName[effectsToInclude]
			dimnames(G) <- dimnames(SF)
			if (!(sienaFitObject$maxlike || sienaFitObject$FinDiff.method))
			{
				D <- derivativeFromScoresAndDeviations(SF, G)
			}
			else
			{
				DF <- sienaFitObject$
						sdf2[ , , effectsToInclude, effectsToInclude,
						drop=FALSE]
				D <- t(apply(DF, c(3, 4), mean))
			}
			fra <- apply(G, 3, sum) / nSims
			doTests <- rep(FALSE, sum(effectsToInclude))
			names(doTests) <- effectsObject$effectName[effectsToInclude]
			doTests[effectsObject$effectName[index]] <- TRUE
			mmThetaDelta <- as.numeric(ScoreTest(length(doTests), D,
							sigma, fra, doTests,
							maxlike=sienaFitObject$maxlike)$oneStep )
			mmPartialThetaDelta <- rep(0,length(theta0))
			mmPartialThetaDelta[length(theta0)] <- mmThetaDelta[length(theta0)]
			JacobianExpStat <- lapply(period, function (i) {
				t(SF[,i,]) %*% simStatsByPeriod[[i]]/ nSims	 })
			Gradient <- lapply(period, function(i) {
						-2	* JacobianExpStat[[i]] %*%
								covInvByPeriod[[i]] %*%
								t( obsStatsByPeriod[[i]] - ExpStat[[i]] ) })
			Hessian <- lapply(period, function (i) {
							2 *
							JacobianExpStat[[i]] %*%
							covInvByPeriod[[i]] %*%
							t(JacobianExpStat[[i]])
					})
			gmmThetaDelta <- -1 * as.numeric( ginv(Reduce("+", Hessian)) %*%
							Reduce("+", Gradient) )
			OneStepSpecs[effectsToInclude,counterTestEffects] <- theta0 +
					mmThetaDelta
			PartialOneStepSpecs[effectsToInclude,counterTestEffects] <-
					theta0 + mmPartialThetaDelta
			GmmOneStepSpecs[effectsToInclude,counterTestEffects] <- theta0 +
					gmmThetaDelta
			for (i in 1:length(obsMhd)) {
				OneStepMHD[[i]][counterTestEffects] <-	as.numeric(
					obsMhd[i] +
					mmThetaDelta %*% Gradient[[i]] + 0.5 *
					mmThetaDelta %*% Hessian[[i]] %*% mmThetaDelta)
				GmmMhdValue[[i]][counterTestEffects] <-
						as.numeric( obsMhd[i] +
						gmmThetaDelta %*%
						Gradient[[i]] + 0.5 *
						gmmThetaDelta %*%
						Hessian[[i]] %*%
						gmmThetaDelta )
				PartialOneStepMHD[[i]][counterTestEffects] <-
						as.numeric(
						obsMhd[i] +
						mmPartialThetaDelta %*%
						Gradient[[i]] +
						0.5 *
						mmPartialThetaDelta %*%
						Hessian[[i]] %*%
						mmPartialThetaDelta)
			}
			JoinedOneStepMHD[counterTestEffects] <-
					Reduce("+",OneStepMHD)[counterTestEffects]
			JoinedPartialOneStepMHD[counterTestEffects] <-
					Reduce("+",PartialOneStepMHD)[counterTestEffects]
			JoinedGmmMhdValue[counterTestEffects] <-
					Reduce("+",GmmMhdValue)[counterTestEffects]
		}
	}

	names(res) <- names(obsStats)
	class(res) <- "sienaGOF"
	attr(res, "scoreTest") <- (sum(sienaFitObject$test) > 0)
	attr(res, "originalMahalanobisDistances") <- obsMhd
	attr(res, "joinedOneStepMahalanobisDistances") <-
			JoinedOneStepMHD
	attr(res, "oneStepSpecs") <- OneStepSpecs
	attr(res, "partialOneStepMahalanobisDistances") <- PartialOneStepMHD
	attr(res, "joinedPartialOneStepMahalanobisDistances") <-
			JoinedPartialOneStepMHD
	attr(res, "partialOneStepSpecs") <- PartialOneStepSpecs
	attr(res, "gmmOneStepSpecs") <- GmmOneStepSpecs
	attr(res, "gmmOneStepMahalanobisDistances") <- GmmMhdValue
	attr(res, "joinedGmmOneStepMahalanobisDistances") <- JoinedGmmMhdValue
	attr(res,"auxiliaryStatisticName") <-
			attr(obsStats,"auxiliaryStatisticName")
	attr(res, "simTime") <- attr(simStats,"time")
	attr(res, "twoTailed") <- twoTailed
	attr(res, "joined") <- join
	res
}

##@print.sienaGOF siena07 Print method for sienaGOF
print.sienaGOF <- function (x, ...) {
	## require(Matrix)
	levls <- 1:length(x)
	pVals <- sapply(levls, function(i) x[[i]]$p)
	titleStr <- "Monte Carlo Mahalanobis distance test p-value: "

	if (! attr(x,"joined"))
	{
		cat("Siena Goodness of Fit (",
			attr(x,"auxiliaryStatisticName"),"),", length(levls)," periods\n=====\n")
		cat(" >",titleStr, "\n")
		for (i in 1:length(pVals))
		{
			cat(names(x)[i], ": ", pVals[i], "\n")
		}
		for (i in 1:length(pVals))
		{
			if (x[[i]]$Rank < dim(x[[i]]$Observations)[2])
			{
				cat(" * Note for", names(x)[i],
					": Only", x[[i]]$Rank, "statistics are",
					"necessary in the auxiliary function.\n")
			}
		}
	}
	else
	{
		cat("Siena Goodness of Fit (",
			attr(x,"auxiliaryStatisticName"),"), all periods\n=====\n")
		cat(titleStr,pVals[1], "\n")
		if (x[[1]]$Rank < dim(x[[1]]$Observations)[2])
			{
				cat("**Note: Only", x[[1]]$Rank, "statistics are",
				"necessary in the auxiliary function.\n")
			}
	}

	if ( attr(x, "twoTailed") )
	{
		cat("-----\nTwo tailed test used.")
	}
	else
	{
		cat("-----\nOne tailed test used ",
		"(i.e. estimated probability of greater distance than observation).\n")
	}
	originalMhd <- attr(x, "originalMahalanobisDistances")
	if (attr(x, "joined")) {
		cat("-----\nCalculated joint MHD = (",
				sum(originalMhd),") for current model.\n")
	}
	else
	{
		for (j in 1:length(originalMhd)) {
			cat("-----\nCalculated period ", j, " MHD = (",
					originalMhd[j],") for current model.\n")
		}
	}
}

##@summary.sienaGOF siena07 Summary method for sienaGOF
summary.sienaGOF <- function(object, ...) {
	x <- object
	print(x)
	if (attr(x, "scoreTest")) {
		oneStepSpecs <- attr(x, "oneStepSpecs")
		oneStepMhd <- attr(x, "oneStepMahalanobisDistances")
# Tom deleted the following statements, because they lead to warnings
# when checking. In the current code they are superfluous.
#		gmmMhd <- attr(x, "gmmOneStepMahalanobisDistances")
#		gmmOneStepSpecs <- attr(x, "gmmOneStepSpecs")
#		partialOneStepSpecs <- attr(x, "partialOneStepSpecs")
#		partialOneStepMhd <- attr(x, "partialOneStepMahalanobisDistances")
#		joinedPartialOneStepMhd <-
#				attr(x, "joinedPartialOneStepMahalanobisDistances")
		joinedOneStepMhd <- attr(x, "joinedOneStepMahalanobisDistances")
#		joinedGmmMhd <- attr(x, "joinedGmmOneStepMahalanobisDistances")
		cat("\nOne-step estimates for modified models.\n")
		if (attr(x, "joined")) {
			for (i in 1:ncol(oneStepSpecs)) {
				a <- cbind(oneStepSpecs[,i, drop=FALSE] )
				b <- matrix( c(joinedOneStepMhd[i] ), ncol=1)
				rownames(b) <- c("MHD")
				a <- rbind(a, b)
				a <- round(a, 3)
				cat("\n**Model", colnames(a)[1], "\n")
				colnames(a) <- "MMD"
				print(a)
			}
		}
		else
		{
			for (j in 1:length(oneStepMhd)) {
				for (i in 1:ncol(oneStepSpecs)) {
					a <- cbind( oneStepSpecs[,i, drop=FALSE] )
					b <- matrix( c(oneStepMhd[[j]][i], ncol=1) )
					rownames(b) <- c("MHD")
					a <- rbind(a, b)
					a <- round(a, 3)
					cat("\n**Model", colnames(a)[1], "\n")
					colnames(a) <- c("MMD")
					print(a)
				}
			}
		}
		cat("\n-----")
	}
	cat("\nComputation time for auxiliary statistic calculations on simulations: ",
			attr(x, "simTime")["elapsed"] , "seconds.\n")
}

##@plot.sienaGOF siena07 Plot method for sienaGOF
plot.sienaGOF <- function (x, center=FALSE, scale=FALSE, violin=TRUE,
		key=NULL, perc=.05, period=1, main=main, ylab=ylab,	 ...)
{
	require(lattice)
	args <- list(...)
	if (is.null(args$main))
	{
		main=paste("Goodness of Fit of",
				attr(x,"auxiliaryStatisticName"))
		if (!attr(x,"joined"))
		{
			main = paste(main, "Period", period)
		}
	}
	else
	{
		main=args
	}

	if (attr(x,"joined"))
	{
		x <- x[[1]]
	}
	else
	{
		x <- x[[period]]
	}
	sims <- x$Simulations
	obs <- x$Observations
	itns <- nrow(sims)
#	vars <- ncol(sims)
	## Need to check for useless statistics here:
	n.obs <- nrow(obs)

	if (center)
	{
		sims.median <- apply(sims, 2, median)
		sims <- sapply(1:ncol(sims), function(i)
					(sims[,i] - sims.median[i]) )
		obs <- matrix(sapply(1:ncol(sims), function(i)
							(obs[,i] - sims.median[i])), nrow=n.obs )
	}
	if (scale)
	{
		sims.min <- apply(sims, 2, min)
		sims.max <- apply(sims, 2, max)
		sims <- sapply(1:ncol(sims), function(i) sims[,i]/(sims.max[i] -
								sims.min[i] ) )
		obs <- matrix(sapply(1:ncol(sims), function(i) obs[,i] /(sims.max[i] -
										sims.min[i] )
				), nrow=n.obs )
	}

	if (is.null(args$ylab))
	{
		ylabel = "Statistic"
		if (center && scale) {
			ylabel = "Statistic (centered and scaled)"
		}
		else if (scale)
		{
			ylabel = "Statistic (scaled)"
		}
		else if (center)
		{
			ylabel = "Statistic (center)"
		}
		else
		{
			ylabel = "Statistic"
		}
	}
	else
	{
		ylabel = args$ylab
	}

	screen <- sapply(1:ncol(obs),function(i){
						(sum(is.nan(rbind(sims,obs)[,i])) == 0) }) &
			(diag(var(rbind(sims,obs)))!=0)
	sims <- sims[,screen, drop=FALSE]
	obs <- obs[,screen, drop=FALSE]
	obsLabels <- round(x$Observations[,screen, drop=FALSE],3)

	if (is.null(args$xlab))
	{
		xlabel = paste( paste("p:", round(x$p, 3),
						collapse = " "), collapse = "\n")
	}
	else
	{
		xlabel = args$xlab
	}

	xAxis <- (1:sum(screen))

	if (is.null(key))
	{
		if (is.null(attr(x, "key")))
		{
			key=xAxis
		}
		else
		{
			key <- attr(x,"key")
		}
	}
	else
	{
		if (length(key) != ncol(obs))
		{
			stop("Key length does not match the number of variates.")
		}
	}

	br <- trellis.par.get("box.rectangle")
	br$col <- 1
	trellis.par.set("box.rectangle", br)
	bu <- trellis.par.get("box.umbrella")
	bu$col <- 1
	trellis.par.set("box.umbrella", bu)
	plot.symbol <- trellis.par.get("plot.symbol")
	plot.symbol$col <- "black"
	plot.symbol$pch <- 4
	plot.symbol$cex <- 1
	trellis.par.set("plot.symbol", plot.symbol)

	panelFunction <- function(..., x=x, y=y, box.ratio){
		ind.lower = max( round(itns * perc/2), 1)
		ind.upper = round(itns * (1-perc/2))
		yperc.lower = sapply(1:ncol(sims), function(i)
					sort(sims[,i])[ind.lower]  )
		yperc.upper = sapply(1:ncol(sims), function(i)
					sort(sims[,i])[ind.upper]  )
		if (violin) {
			panel.violin(x, y, box.ratio=box.ratio, col = "transparent", ...)
		}
		panel.bwplot(x, y, box.ratio=.1, fill = "gray", ...)
		panel.xyplot(xAxis, yperc.lower, lty=3, col = "gray", lwd=3, type="l",
				...)
		panel.xyplot(xAxis, yperc.upper, lty=3, col = "gray", lwd=3, type="l",
				...)
		for(i in 1:nrow(obs))
		{
			panel.xyplot(xAxis, obs[i,],  col="red", type="l", lwd=1, ...)
			panel.xyplot(xAxis, obs[i,],  col="red", type="p", lwd=3, pch=19,
					...)
			panel.text(xAxis, obs[i,], labels=obsLabels[i,], pos=4)
		}
	}
	bwplot(as.numeric(sims)~rep(xAxis, each=itns), horizontal=FALSE,
			panel = panelFunction, xlab=xlabel, ylab=ylabel,
			scales=list(x=list(labels=key), y=list(draw=FALSE)),
			main=main, ...)
}

##@sparseMatrixExtraction sienaGOF Extracts simulated networks
# This function returns the simulated network as a dgCMatrix;
# this is the "standard" class for sparse numeric matrices
# in the Matrix package. See the help file for "dgCMatrix-class".
# Ties for ordered pairs with a missing value for wave=period or period+1
#  are zeroed;
# note that this also is done in RSiena for calculation of target statistics.
sparseMatrixExtraction <-
	function(i, obsData, sims, period, groupName, varName){
	# require(Matrix)
	dimsOfDepVar<- attr(obsData[[groupName]]$depvars[[varName]], "netdims")
	if (attr(obsData[[groupName]]$depvars[[varName]], "sparse"))
	{
		missings <-
			(is.na(obsData[[groupName]]$depvars[[varName]][[period]]) |
			is.na(obsData[[groupName]]$depvars[[varName]][[period+1]]))*1
	}
	else
	{
		missings <- Matrix(
			(is.na(obsData[[groupName]]$depvars[[varName]][,,period]) |
			is.na(obsData[[groupName]]$depvars[[varName]][,,period+1]))*1)
	}
	if (is.null(i))
	{
		# sienaGOF wants the observation;
		# transform structurally fixed values into regular values
		# by "modulo 10" operation
		if (attr(obsData[[groupName]]$depvars[[varName]], "sparse"))
		{
			returnValue <- drop0(
			 Matrix(obsData[[groupName]]$depvars[[varName]][[period+1]] %% 10))
		}
		else
		{
			returnValue <-
			 Matrix(obsData[[groupName]]$depvars[[varName]][,,period+1] %% 10)
		}
		returnValue[is.na(returnValue)] <- 0
	}
	else
	{
		# sienaGOF wants the i-th simulation:
		returnValue <- sparseMatrix(
				sims[[i]][[groupName]][[varName]][[period]][,1],
				sims[[i]][[groupName]][[varName]][[period]][,2],
				x=sims[[i]][[groupName]][[varName]][[period]][,3],
				dims=dimsOfDepVar[1:2] )
	}
	## Zero missings:
	1*((returnValue - missings) > 0)
}

##@networkExtraction sienaGOF Extracts simulated networks
# This function provides a standard way of extracting simulated and observed
# networks from the results of a siena07 run.
# It returns the network as an edge list of class "network"
# according to the <network> package (used for package sna).
# Ties for ordered pairs with a missing value for wave=period or period+1
# are zeroed;
# note that this also is done in RSiena for calculation of target statistics.
networkExtraction <- function (i, obsData, sims, period, groupName, varName){
	require(network)
	dimsOfDepVar<- attr(obsData[[groupName]]$depvars[[varName]], "netdims")
	isbipartite <- (attr(obsData[[groupName]]$depvars[[varName]], "type")
						=="bipartite")
	sparseData <- (attr(obsData[[groupName]]$depvars[[varName]], "sparse"))
	# For bipartite networks in package <network>,
	# the number of nodes is equal to
	# the number of actors (rows) plus the number of events (columns)
	# with all actors preceeding all events.
	# Therefore the bipartiteOffset will come in handy:
	bipartiteOffset <- ifelse (isbipartite, 1 + dimsOfDepVar[1], 1)

	# Initialize empty networks:
	if (isbipartite)
	{
		# for bipartite networks, package <<network>> numbers
		# the second mode vertices consecutively to the first mode.
		emptyNetwork <- network.initialize(dimsOfDepVar[1]+dimsOfDepVar[2],
											bipartite=dimsOfDepVar[1])
	}
	else
	{
		emptyNetwork <- network.initialize(dimsOfDepVar[1],	bipartite=NULL)
	}
	if (sparseData)
	{
		# Which tie variables are regarded as missings:
		missings <- as(
			is.na(obsData[[groupName]]$depvars[[varName]][[period]]) |
			is.na(obsData[[groupName]]$depvars[[varName]][[period+1]]),
			"lgTMatrix")
		# For lgTMatrix, slots i and j are the rows and columns,
		# numbered from 0 to dimension - 1.
		# Actors in class network are numbered starting from 1.
		# Hence 1 must be added to missings@i and missings@j.
		# Put the missings into network shape:
		if (length(missings@i) <= 0)
		{
			missings <- emptyNetwork
		}
		else
		{
			missings <- network.edgelist(
			  cbind(missings@i + 1, missings@j + bipartiteOffset, 1), emptyNetwork)
		}
	}
	else # not sparse
	{
		# For adjacency matrices the size is evident.
		if (isbipartite)
		{
			missings <- network(
				(is.na(obsData[[groupName]]$depvars[[varName]][,,period]) |
				is.na(obsData[[groupName]]$depvars[[varName]][,,period+1]))*1,
				matrix.type="adjacency", bipartite=dimsOfDepVar[1])
		}
		else
		{
			missings <- network(
				(is.na(obsData[[groupName]]$depvars[[varName]][,,period]) |
				is.na(obsData[[groupName]]$depvars[[varName]][,,period+1]))*1,
				matrix.type="adjacency")
		}
	}

	if (is.null(i))
	{
		# sienaGOF wants the observation;
		if (sparseData)
		{
		# transform structurally fixed values into regular values
		# by "modulo 10" operation;
		# drop NAs and 0 values
			original <-
				obsData[[groupName]]$depvars[[varName]][[period+1]] %% 10
			original[is.na(original)] <- 0
			original <- as(drop0(original), "dgTMatrix")
		# now original@x is a column of ones;
		# the 1 here is redundant because of the default ignore.eval=TRUE
		# in network.edgelist
			returnValue <- network.edgelist(
					cbind(original@i + 1, original@j + bipartiteOffset, 1),
					emptyNetwork) - missings

		}
		else # not sparse: deal with adjacency matrices
		{
			original <-
				obsData[[groupName]]$depvars[[varName]][,,period+1] %% 10
			original[is.na(original)] <- 0
			if (isbipartite)
			{
				returnValue <- network(original, matrix.type="adjacency",
									bipartite=dimsOfDepVar[1]) - missings
			}
			else
			{
				returnValue <- network(original, matrix.type="adjacency",
									bipartite=FALSE) -	missings
			}
		}
	}
	else
	{
		# sienaGOF wants the i-th simulation;
		# the 1 in cbind is redundant because of the default ignore.eval=TRUE
		# in network.edgelist
		bipartiteOffset <- ifelse (isbipartite, dimsOfDepVar[1], 0)
		returnValue <- network.edgelist(
			cbind(
			 sims[[i]][[groupName]][[varName]][[period]][,1],
			(sims[[i]][[groupName]][[varName]][[period]][,2] + bipartiteOffset),
						1), emptyNetwork) - missings
	}
  returnValue
}

##@behaviorExtraction sienaGOF Extracts simulated behavioral variables.
# This function provides a standard way of extracting simulated and observed
# dependent behavior variables from the results of a siena07 run.
# The result is an integer vector.
# Values for actors with a missing value for wave=period or period+1 are
# transformed to NA.
behaviorExtraction <- function (i, obsData, sims, period, groupName, varName) {
  missings <- is.na(obsData[[groupName]]$depvars[[varName]][,,period]) |
	is.na(obsData[[groupName]]$depvars[[varName]][,,period+1])
  if (is.null(i))
	{
		# sienaGOF wants the observation:
		original <- obsData[[groupName]]$depvars[[varName]][,,period+1]
		original[missings] <- NA
		returnValue <- original
	}
	else
	{
		#sienaGOF wants the i-th simulation:
		returnValue <- sims[[i]][[groupName]][[varName]][[period]]
		returnValue[missings] <- NA
	}
	returnValue
}

##@OutdegreeDistribution sienaGOF Calculates Outdegree distribution
OutdegreeDistribution <- function(i, obsData, sims, period, groupName, varName,
						levls=0:8, cumulative=TRUE) {
	x <- sparseMatrixExtraction(i, obsData, sims, period, groupName, varName)
	a <- apply(x, 1, sum)
	if (cumulative)
	{
		oddi <- sapply(levls, function(i){ sum(a<=i) })
	}
	else
	{
		oddi <- sapply(levls, function(i){ sum(a==i) })
	}
	names(oddi) <- as.character(levls)
	oddi
}

##@IndegreeDistribution sienaGOF Calculates Indegree distribution
IndegreeDistribution <- function (i, obsData, sims, period, groupName, varName,
						levls=0:8, cumulative=TRUE){
  x <- sparseMatrixExtraction(i, obsData, sims, period, groupName, varName)
  a <- apply(x, 2, sum)
  if (cumulative)
  {
	iddi <- sapply(levls, function(i){ sum(a<=i) })
  }
  else
  {
	iddi <- sapply(levls, function(i){ sum(a==i) })
  }
  names(iddi) <- as.character(levls)
  iddi
}

##@BehaviorDistribution sienaGOF Calculates behavior distribution
BehaviorDistribution <- function (i, obsData, sims, period, groupName, varName,
							levls=NULL, cumulative=TRUE){
	x <- behaviorExtraction(i, obsData, sims, period, groupName, varName)
	if (is.null(levls))
	{
		levls <- attr(obsData[[groupName]]$depvars[[varName]],"behRange")[1]:
					attr(obsData[[groupName]]$depvars[[varName]],"behRange")[2]
	}
	if (cumulative)
	{
		bdi <- sapply(levls, function(i){ sum(x<=i, na.rm=TRUE) })
	}
	else
	{
	bdi <- sapply(levls, function(i){ sum(x==i, na.rm=TRUE) })
	}
	names(bdi) <- as.character(levls)
	bdi
}