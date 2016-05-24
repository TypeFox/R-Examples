impulseResponseVAR1 <- function(A, Tmax=10, figure=FALSE){
	######################################################################################
	#
	# DESCRIPTION:
	# Calculate impulse responses for VAR(1) model. MA representation of the VAR(1) 
	# model dictates that these are given by powers of A.
	#
	# ARGUMENTS:
	# -> A       : Matrix A of regression parameters.
	# -> Tmax    : Maximum mumber of time points up to which the impulse responses are to be evaluated.
	# -> figure  : Logical, indicating whether a summary plot of the impulse responses should be generated.
	#
	# DEPENDENCIES:
	# library(expm)	    # functions: %^% 
	#
	# NOTES:
	# ...
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

	# calculate impulse responses
	impResponse <- apply(abs(A), 2, mean)
	for (t in 2:Tmax){
		impResponse <- cbind(impResponse, apply(abs(A %^% t), 2, mean))
	}
    
	# assess top 5: covariates with largest (absolute) impulse response
	ranks <- rank(rowSums(impResponse))
	top <- which(ranks > nrow(A) - min(nrow(A), 5))
	ranks <- nrow(A) - ranks + 1
    
	# plot absolute impulse response against time
	if (figure){
		cols <- c("red4", "red", "orangered", "darkorange1", "orange")
		plot(impResponse[1,] ~ c(1:Tmax), xlab="time", ylab="mean absolute impulse response", ylim=c(0, max(impResponse)), col="white", type="l", lty=1, main="impulse response analysis")
		for (j in 1:nrow(A)){
			if (j %in% top){
				lines(impResponse[j,] ~ c(1:Tmax), col=cols[ranks[j]], lty=ranks[j], lwd=2)
			} else {
				lines(impResponse[j,] ~ c(1:Tmax), col="grey", lty=ranks[j], lwd=1)
			}
        	}
	    	legend("topright", paste("variate ", as.character(top), sep=""), col=cols[ranks[top]], lty=ranks[top], lwd=2)		
    	}
    
	return(impResponse)   
}

