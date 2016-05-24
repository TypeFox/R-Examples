cdfDist <-
function(x1, F1, x2, F2){
	if(length(unique(x1)) != length(x1)){
		stop("The values in x1 are not unique")
	} else {
		temp <- order(x1)
		obsX <- x1[temp]
		F1 <- F1[temp]
	}
	if(length(unique(x2)) != length(x2)){
		stop("The values in x2 are not unique")
	} else {
		temp <- order(x2)
		expX <- x2[temp]
		F2 <- F2[temp]
	}
	mu    <- 0
	xi    <- sort(union(obsX, expX))
	temp  <- max(c(obsX[1], expX[1]))
	temp1 <- min(c(tail(obsX, 1), tail(expX, 1)))
	xi    <- xi[xi >= temp & xi <= temp1]
	n     <- length(xi)
	F1    <- F1[obsX >= temp & obsX <= temp1]
	obsX  <- obsX[obsX >= temp & obsX <= temp1]
	F2    <- F2[expX >= temp & expX <= temp1]
	expX  <- expX[expX >= temp & expX <= temp1]
  #	xiObs <- which(obsQ %in% xi)
  #	xiExp <- which(expQ %in% xi)
	if(sum(obsX %in% xi) == 0){
		stop("")
	} else {
		oF               <- rep(NA, n)
		oF[xi %in% obsX] <- F1[obsX %in% xi]
	}
	if(sum(expX %in% xi) == 0){
		stop("")
	} else {
		eF               <- rep(NA, n)
		eF[xi %in% expX] <- F2[expX %in% xi]
	}
	for(i in which(is.na(oF))){
		oF[i] <- F1[which.min(abs(obsX - xi[i]))]
	}
	for(i in which(is.na(eF))){
		eF[i] <- F2[which.min(abs(expX - xi[i]))]
	}
	d    <- cbind(xi, oF, eF)
	maxF <- apply(d[,2:3], 1, max)
	mu   <- diff(xi)*(oF[-1]-eF[-1])^2/(1-maxF[-1])/diff(range(xi))
	mu   <- cumsum(mu)
	xi   <- xi[-1]
	tR      <- list(x=xi,
					F1=oF[-1],
					F2=eF[-1],
					meas=mu,
					cdfDist=tail(mu,1))
	class(tR) <- "cdfDist"
	return(tR)
}

