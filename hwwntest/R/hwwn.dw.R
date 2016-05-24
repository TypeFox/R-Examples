hwwn.dw <-
function(J, filter.number, filter.family) {

# Used by nf.thresholds
	
	longest.support <- (2^J - 1)*(2*filter.number - 1) + 1
	J1 <- ceiling(logb(longest.support, 2))
	tmp <- rep(0, max(2^J1,4))
	tmp.w <- wd(tmp, filter.number, filter.family)
	dw <- vector("list",  J)
	index <- filter.number
	for (j in 1:J) {
		tmp.w$D[index] <- 1
		dw[[j]] <- wr(tmp.w)[1:((2^j - 1)*(2*filter.number - 1) + 1)]
		tmp.w$D[index] <- 0
		index <- index + 2^(J1-j)
	}
	return(dw)
}

