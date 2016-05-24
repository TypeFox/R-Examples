"escouf" <-
function(x, level=1, verbose=TRUE) {
	call <- match.call()
	# Rem: we could decide to store the initial data into res$data
	# To free memory, we will just store a call to these data
	# The drawback is that initial data should not be modified
	# between 'ecouf' and 'extract.escouf'!!!
	Data <- deparse(substitute(x))
	# Calculate the trace of a matrix (sum of its diagonal elements)
	Trace <- function(x) {sum(c(x)[1 + 0:(min(dim(x)) - 1) * (dim(x)[1] + 1)], na.rm=TRUE)}
		
	x <- as.data.frame(x)	# We want to be sure to work on a data frame!
	Names <- names(x)
	p <- ncol(x)
	vt <- 1:p					# Variable to test
	vr <- NULL					# Held variables
	vrt <- NULL					# Temporarily held variables
	RV <- NULL					# Final held variables
	for (i in 1:p) {			# Loop on the number of variables
		Rvmax <- 0
		for (j in 1:(p-i+1)) {	# loop on variables
			if (!is.null(vr)) {	# New table
				x2 <- cbind(x, x[vr], x[vt[j]])
			} else {
				x2 <- cbind(x, x[vt[j]])
			}	
			Rtot <- cor(x2)		# Correlations table
			Ryy <- Rtot[1:p, 1:p]
			Rxx <- Rtot[(p+1):(p+i), (p+1):(p+i)]
			Rxy <- Rtot [(p+1):(p+i), 1:p]
			Ryx <- t(Rxy)
			rv <- Trace(Ryx %*% Rxy)/sqrt(Trace(Ryy %*% Ryy)*Trace(Rxx %*% Rxx))	# rv calculation
			if (rv>Rvmax) {
				Rvmax <- rv		# Test on rv
				vrt <- vt[j]	# Temporarily held variable
			}
		}
		vr[i] <- vrt
		vt <- vt[vt!=vr[i]]		# Reidentify variables to test
		RV[i] <- Rvmax			# Final held variable
		if (verbose==TRUE) {
			vrStr <- format(c(vr[i], 111))[1]
			cat("Variable", vrStr, "incorporated, RV =", Rvmax, "\n")
			flush.console()
		}
		if (Rvmax>level) break	# Stop iteration (level reached)
	}
	names(vr) <- Names[vr]		# Gives variable names to vr
	names(RV) <- Names[vr]		# ... and to RV
	res <- list(data=Data, vr=vr, RV=RV, calc.level=level, vars=c(p, length(vr)), call=call)		# Create a list containing the result
	class(res) <- "escouf"		# and turn it into an 'escouf' object
	res							# Return the result
}
