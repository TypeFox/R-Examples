mutualInfoVAR1 <- function(A, SigmaE, Tmax=10, figure=FALSE){
	######################################################################################
	#
	# DESCRIPTION:
	# Calculate the mutual informations for VAR(1) model.
	#
	# ARGUMENTS:
	# -> A        : Matrix A of regression parameters.
	# -> SigmaE   : Covariance matrix of the errors (innovations).
	# -> Tmax     : Maximum mumber of time points up to which the mutual informations are to be evaluated.
	# -> figure   : Logical, indicating whether a summary plot of the mutual informations should be generated.    
	# 
	# DEPENDENCIES:
	# library(expm)	    # functions: %^%
	#
	# NOTES:
	# ....
	# 	
	######################################################################################

	# input checks
	if (as.character(class(A)) != "matrix"){ stop("Input (A) is of wrong class.") }
	if (nrow(A) != ncol(A)){ stop("Matrix A is not square.") }
	if (as.character(class(Tmax)) != "numeric"){ stop("Input (Tmax) is of wrong class.") }
	if (length(Tmax) != 1){ stop("Input (Tmax) is of wrong length.") }
	if (is.na(Tmax)){ stop("Input (Tmax) is not a positive integer.") }
	if (Tmax < 1){ stop("Input (Tmax) should be an integer larger than one.") }
	if (as.character(class(figure)) != "logical"){ stop("Input (figure) is of wrong class.") }
	if (as.character(class(SigmaE)) != "matrix"){ stop("Input (SigmaE) is of wrong class.") }    
	if (!isSymmetric(SigmaE)){ stop("Non-symmetrical matrix for the error covariance matrix provided") } 
	if (nrow(A) != nrow(SigmaE)){ stop("Dimensions of input (A) do not match that of other input (SigmaE).") } 
	if (ncol(A) != ncol(SigmaE)){ stop("Dimensions of input (A) do not match that of other input (SigmaE).") } 
    
	# calculate impulse responses
	MI <- function(j, A, SigmaE, Tmax){
		# covariance matrix with right zero
		SigmaEnull <- SigmaE
		SigmaEnull[j,] <- SigmaEnull[,j] <- 0
		SigmaEnull[-j,-j] <- SigmaEnull[-j,-j] - SigmaE[-j, j, drop=FALSE] %*% solve(SigmaE[j,j]) %*% SigmaE[j,-j,drop=FALSE]

		# initial value of marginal and conditional covariances
		MIslh <- numeric()
		varMarg <- varCond <- SigmaE 
		for (tau in 1:Tmax){
			Atau <- A %^% tau
			varCond <- varMarg + Atau %*% SigmaEnull %*% t(Atau)
			varMarg <- varMarg + Atau %*% SigmaE %*% t(Atau)
 			MIslh <- c(MIslh, log(det(varMarg)) - log(det(varCond)))
		}
		return(MIslh)
	}    
	MIs <- matrix(unlist(lapply(1:nrow(A), MI, A=A, SigmaE=SigmaE, Tmax=Tmax)), ncol=Tmax, byrow=TRUE)
         
	# assess top 5: covariates with largest mutual informations
	ranks <- rank(rowSums(MIs))
	top <- which(ranks > nrow(A) - min(nrow(A), 5))
	ranks <- nrow(A) - ranks + 1
    
	# plot absolute impulse response against time
	if (figure){
		cols <- c("red4", "red", "orangered", "darkorange1", "orange")
		plot(MIs[1,] ~ c(1:Tmax), xlab="time", ylab="mutual information", ylim=c(0, max(MIs)), col="white", type="l", lty=1, main="mutual information")
		for (j in 1:nrow(A)){
			if (j %in% top){
				lines(MIs[j,] ~ c(1:Tmax), col=cols[ranks[j]], lty=ranks[j], lwd=2)
			} else {
				lines(MIs[j,] ~ c(1:Tmax), col="grey", lty=ranks[j], lwd=1)
			}
        	}
	    	legend("topright", paste("variate ", as.character(top), sep=""), col=cols[ranks[top]], lty=ranks[top], lwd=2)		
    	}
        
    return(MIs)   
}


