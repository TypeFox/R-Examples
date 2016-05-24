ridgePathVAR1 <- function(Y, lambdaAgrid, lambdaPgrid, pathType="A", plotTypeSigmaE="pcor", diag=FALSE, verbose=TRUE, ...){
	#####################################################################################################
	# 
	# DESCRIPTION:	
	# Function that visualizes the regularization paths of the parameters of the VAR(1) model. The elements of
	# the ridge ML estimate of either A or SigmaE are plotted against a specified range of their penalty 
	# parameter (keeping the other penalty parameter fixed). 
	#
	# 
	# ARGUMENTS:
	# -> Y              : Three-dimensional array containing the data. The first, second and third dimensions correspond to 
	#                     covariates, time and samples, respectively. The data are assumed to centered covariate-wise.
	# -> lambdaAgrid    : Numeric of length larger than one. It contains the grid points corresponding to the lambdaA.
	# -> lambdaPgrid    : Numeric of length larger than one. It contains the grid points corresponding to the lambdaO.
	# -> pathType       : character indicating of which parameter to plot its ridge regularization paths. Either "A" or "SigmaE". 
	# -> plotTypeSigmaE : A character indicating the type of element for which a visualization of the regularization paths 
	#                     (of SigmaE) is desired. Must be one of: "pcor", "cor", "cov", "prec". 
	# -> diag           : A logical indicating if the diagonal elements should be retained for visualization of the 
	#                     regularization path of SigmaE.
	# -> verbose        : Logical indicator: should intermediate output be printed on the screen?
	# -> ...            :  Other arguments to be passed to \code{ridgeVAR1}.
	# 
	# DEPENDENCIES:
	# library(rags2ridges)	    # functions: ridgeS, ridgeSchordal, default.target
	#
	# NOTES: 
	# .... 
	#
	################################################################

	# input checks
	if (class(verbose) != "logical"){ stop("Input (verbose) is of wrong class") }
	if (verbose){ cat("Perform input checks...", "\n")   }
	if (as.character(class(Y)) != "array"){ stop("Input (Y) is of wrong class.") }
	if (length(dim(Y)) != 3){ stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.") }
	if (class(lambdaAgrid) != "numeric"){ stop("Input (lambdaAmin) is of wrong class") }
	if (length(lambdaAgrid) < 2){ stop("lambdaAgrid of length two or longer") }
	if (any(lambdaAgrid <= 0)){ stop("all values of lambdaAgrid must be positive") }
	if (class(lambdaPgrid) != "numeric"){ stop("Input (lambdaAmin) is of wrong class") }
	if (length(lambdaPgrid) < 2){ stop("lambdaAgrid of length two or longer") }
	if (any(lambdaPgrid <= 0)){ stop("all values of lambdaAgrid must be positive") }
	if (class(pathType) != "character"){ stop("Input (pathType) is of wrong class") }
	if (length(intersect(pathType, c("A", "SigmaE"))) == 0){ stop("Input (pathType) is of wrongly specified") }
	if (class(plotTypeSigmaE) != "character"){ stop("Input (plotTypeSigmaE) is of wrong class") }
	if (class(diag) != "logical"){ stop("Input (diag) is of wrong class") }
	if (length(intersect(plotTypeSigmaE, c("pcor", "cor", "cov", "prec"))) == 0){ stop("Input (plotTypeSigmaE) is of wrongly specified") }

	if (verbose){ cat("Calculating...", "\n") }
	if (pathType=="A"){
        lambdaAgrid <- sort(lambdaAgrid)
		lambdaP <- min(lambdaPgrid)
		YforPlot <- numeric()
		for (k in 1:length(lambdaAgrid)) {
			YforPlot <- cbind(YforPlot, as.numeric(ridgeVAR1(Y, lambdaA=lambdaAgrid[k], lambdaP=lambdaP, ...)$A))
			if (verbose) {
				cat(rep("\b", 100), sep="")
				cat(paste("grid point ", k, " (lambdaA = ", lambdaAgrid[k], ") done", sep = ""))
			}
		}
		if (verbose) {	cat("\n") }
		# plot regularization path for A
		plot(YforPlot[1,] ~ log(lambdaAgrid), xlab = expression(log(lambda[a])), ylab = "penalized elements of A", 
			main=bquote(paste("ridge regularization path of A (with ", lambda[omega], "=", .(round(lambdaP, 3)), ")", sep="")), col="white", ylim=c(min(YforPlot), max(YforPlot)))
		for (k in 1:nrow(YforPlot)){ 
			lines(YforPlot[k,] ~ log(lambdaAgrid), col=k, lty=k)
		}
	}

	if (pathType=="SigmaE"){
	        lambdaA <- min(lambdaAgrid)
		lambdaPgrid <- sort(lambdaPgrid)
		YforPlot <- numeric()
		for (k in 1:length(lambdaPgrid)) {
			S <- solve(ridgeVAR1(Y, lambdaA=lambdaA, lambdaP=lambdaPgrid[k])$P)
			if (plotTypeSigmaE=="pcor"){
				YforPlot <- cbind(YforPlot, -cov2cor(solve(S))[upper.tri(S)])
			}
			if (plotTypeSigmaE=="prec"){
				YforPlot <- cbind(YforPlot, solve(S)[upper.tri(S, diag=diag)])
			}
			if (plotTypeSigmaE=="cov"){
				YforPlot <- cbind(YforPlot, S[upper.tri(S, diag=diag)])
			}
			if (plotTypeSigmaE=="cor"){
				YforPlot <- cbind(YforPlot, cov2cor(S)[upper.tri(S)])
			}
			if (verbose) {
				cat(rep("\b", 100), sep="")
				cat(paste("grid point ", k, " (lambdaP = ", lambdaPgrid[k], ") done", sep = ""))
			}
		}
		if (verbose) {	cat("\n") }

		# plot regularization path for SigmaE
	        if (plotTypeSigmaE=="cor"){ ylabel <- "penalized correlation" }
	        if (plotTypeSigmaE=="cov"){ ylabel <- "penalized covariances" }
	        if (plotTypeSigmaE=="pcor"){ ylabel <- "penalized partial correlation" }
	        if (plotTypeSigmaE=="prec"){ ylabel <- "penalized precision elements" }
	        plot(YforPlot[1,] ~ log(lambdaPgrid), xlab = expression(log(lambda[omega])), ylab=ylabel, main=bquote(paste("ridge regularization path of ", Sigma[epsilon], " (with ", lambda[a], "=", .(round(lambdaA, 3)), ")", sep="")), col="white", ylim=c(min(YforPlot), max(YforPlot)))
		for (k in 1:nrow(YforPlot)){ 
			lines(YforPlot[k,] ~ log(lambdaPgrid), col=k, lty=k)
		}
        }
}



