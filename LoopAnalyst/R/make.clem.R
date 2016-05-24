make.clem <- function(CM, status=FALSE) {
	# Make births matrix
	CM.B <- CM
	CM.B[CM.B<0] <- 0

	# Make deaths matrix
	CM.D <- (CM.B - CM)	

	Adj <- make.adjoint(CM, status)

	clem.B <- sign(-CM.B%*%Adj)
	clem.D <- sign(-CM.D%*%Adj)
	
	CLEM <- clem.B

	CLEM[is.na(CLEM)] <- " ? "
	CLEM[CLEM == 1] <- " + "
	CLEM[CLEM == -1] <- " - "
	CLEM[CLEM == 0] <- " 0 "

	clem.D[is.na(clem.D)] <- "?"
	clem.D[clem.D == 1] <- "+"
	clem.D[clem.D == -1] <- "-"
	clem.D[clem.D == 0] <- "0"

	N <- nrow(CLEM)
	for (i in 1:N) {
		if (is.na(clem.B[i,i])) {
			CLEM[i,i] <- "?"
			}
		else if (clem.B[i,i] == 1) {
			CLEM[i,i] <- "+"
			}
		else if (clem.B[i,i] == -1) {
			CLEM[i,i] <- "-"
			}
		else if (clem.B[i,i] == 0) {
			CLEM[i,i] <- "0"
			}
		CLEM[i,i] <- paste(clem.D[i,i],",",CLEM[i,i],sep="")
		}

	dimnames(CLEM)[[2]] <- paste(" ",dimnames(CLEM)[[2]]," ",sep="")

	print(CLEM, quote=FALSE)
	
	}
