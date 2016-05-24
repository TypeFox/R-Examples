du.adjcat <-
function(eta, x, k, Ymat) {
		ceta <- .Call("do_row_cumsum", eta)
	 	eta.cumsum <- matrix(ceta, 
				nrow=nrow(eta),
				byrow=T)
		class1 <- rowSums(col(eta.cumsum) * exp(eta.cumsum)) /
			(1 + rowSums(exp(eta.cumsum)))
		multiplier <- matrix(nrow=dim(x)[1], ncol=k)
		for (j in 1:k) {
			 multiplier[,j] <- (ncol(as.matrix(Ymat[,1:j]))-1) *
			 	Ymat[,j] - Ymat[,j] * class1
		}
		u <- t(x) %*% multiplier
        update.value <- min(rowSums(-u, na.rm=TRUE))
        update.j <- which.min(rowSums(-u, na.rm=TRUE))  
        list(update.j=update.j, update.value=update.value)  
}
