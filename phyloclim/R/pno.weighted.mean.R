pno.weighted.mean <- 						
function(x, subset = NULL, normalize = TRUE){
	
	#x <- cbind("VAR" = as.numeric(rownames(x)), x)
	
	# subset matrix
	# ---------------------
	if (!is.null(subset))
		x <- x[, c(1, which(names(x) %in% subset))]
		
	nb <- dim(x)[2]
	wm <- vector(length = (nb - 1))
	for (i in 2:nb){
		if (normalize){
			nf <- sum(x[, i])
			wm[i - 1] <- sum(x[, 1] * (x[, i] / nf))
		}
		else
			wm[i - 1] <- sum(x[, 1] * x[, i])
	}
	names(wm) <- colnames(x)[2:nb]
	wm
}