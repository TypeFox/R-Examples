impute <- function(RRMatrix, na.rm="meansNA", stress.max = 0.01, maxIt=100) {

	# in Matrix umwandeln, sonst geht's nicht ...
	RRMatrix2 <- as.matrix(RRMatrix)
	rownames(RRMatrix2) <- rownames(RRMatrix)
	colnames(RRMatrix2) <- colnames(RRMatrix)
	RRMatrix <- RRMatrix2
	
	NAs <- is.na(RRMatrix)		# Matrix, die die Position der NAs ausserhalb der Diagonale speichert
	
	save.diag <- diag(RRMatrix)	# self ratings aus der Diagonale abspeichern
	diag(RRMatrix) <- NA

	eff0 <- eff <- quickeffects(RRMatrix)
	stress <- 1
	it <- 0

	# save evolution of parameters
	as <- matrix(eff$a, nrow=1, ncol=ncol(RRMatrix))
	bs <- matrix(eff$b, nrow=1, ncol=ncol(RRMatrix))
	ms <- c()
	
	while (stress > stress.max) {
	
		rM <- rowMeans(RRMatrix, na.rm=TRUE)
		cM <- colMeans(RRMatrix, na.rm=TRUE)
		rM_mean <- mean(rM)
		cM_mean <- mean(cM)
		grandmean <- mean(RRMatrix, na.rm=TRUE)
		
		for (i in 1:ncol(RRMatrix)) {
			for (j in 1:nrow(RRMatrix)) {
				if (NAs[j,i]==TRUE) {
					
					# Ersetzung durch mittleres Zeilen/ Spaltenmittel
					if (grepl("mean", na.rm)) {
						RRMatrix[j,i] <- (rM[j] + cM[i])/2
					}
				
				}
			}
		}

		if (grepl("1", na.rm)) break; # bei means1, chi1: beim ersten Durchgang gleich rausspringen

		eff <- quickeffects(RRMatrix)

		stress <- max(max(abs(eff$a - eff0$a)), max(abs(eff$b - eff0$b)))
		eff0 <- eff
		it <- it + 1
		
		if (it > maxIt) {
			warning("Maximum iterations exceeded; fall back to single imputation.", call.=FALSE)
			return(impute(RRMatrix2, paste(na.rm,"1",sep="")))
		}
		
		as <- rbind(as, eff$a)
		bs <- rbind(bs, eff$b)
		ms <- c(ms, eff$m)
	}	
	
	diag(RRMatrix) <- save.diag
	if (!grepl("NA", na.rm)) {NAs[] <- FALSE}
	
	return(list(RRMatrix=RRMatrix, NAs=NAs, iterations=it, as=as, bs=bs, ms=ms))
}