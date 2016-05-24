`hankel` <-
function(y, lag = 1, cutoff = .90, type = "median") {
    	if(lag != 1 && lag != 2 && lag != 3) {
        	stop("ERROR: ", paste(sQuote("lag"), sep = ""), " must be 1, 2, or 3")
    	}
    	if (type != "median" && type != "mean") {
        	stop("ERROR: ", paste(sQuote("type"), sep = ""), " must be ", 
            	paste(dQuote("median"), sep = ""), " or ", paste(dQuote("mean"), 
                	sep = ""), ".")
    	}
    	if(is.list(y) != TRUE) {
        	stop("Error: ", paste(sQuote("y"), sep = ""), " must be a list.")
    	}
	R <- length(y)
	T <- dim(y[[1]])[2]
	P <- dim(y[[1]])[1]
	svs <- vector("list", R)
	dim <- rep(0, R)
	dim2 <- rep(0,R)
	for(r in 1:R) {
		H <- matrix(0, nrow = lag*P, ncol = lag*P)
		Gamma <- vector("list", 2*lag-1)
		for(i in 1:(2*lag-1)) {
			Gamma[[i]] <- matrix(0, nrow = P, ncol = P)
			for(t in 1:(T-i)) {
				Gamma[[i]] <- Gamma[[i]] + 
					matrix(y[[r]][,t], nrow = P, ncol = 1) %*%
					matrix(y[[r]][,(t+i)], nrow = 1, ncol = P)
			}
		}
		H[1:P, 1:P] <- Gamma[[1]]
		if(lag >= 2) {
			H[(P+1):(2*P), 1:P] <- H[1:P, (P+1):(2*P)] <- Gamma[[2]]
			H[(P+1):(2*P), (P+1):(2*P)] <- Gamma[[3]]
		}
		if(lag == 3) {
			H[(2*P+1):(3*P), 1:(P)] <- H[1:P, (2*P+1):(3*P)] <- Gamma[[3]]
			H[(2*P+1):(3*P), (P+1):(2*P)] <- 
				H[(P+1):(2*P), (2*P+1):(3*P)] <- Gamma[[4]]
			H[(2*P+1):(3*P), (2*P+1):(3*P)] <- Gamma[[5]]
		}
		svs[[r]] <- svd(H)$d/svd(H)$d[1]
		dim[r] <- sumFunc(svs[[r]], cutoff)
	}
	dim.choice <- round(apply(matrix(dim, nrow = 1), 1, type))
	return(list(svs = svs, dim = dim.choice))
}

