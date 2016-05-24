as.data.frame.inflmlm <-
function(x, ..., FUN=det, funnames=TRUE) {
	m <- x$m
	if(m==1) {
		df <- with(x, data.frame(H, Q, CookD, L, R))
		rownames(df) <- x$labels
		}
	else {
#		FUN <- match.arg(FUN)
		H <- sapply(x$H, FUN)
		Q <- sapply(x$Q, FUN)
		L <- sapply(x$L, FUN)
		R <- sapply(x$R, FUN)
		df <- data.frame(H, Q, CookD=x$CookD, L, R)
		rownames(df) <- apply(x$subsets,1, paste, collapse=',')
		if(funnames) colnames(df)[c(1,2,4,5)] <- paste(deparse(substitute(FUN)), c("H","Q","L","R"), sep="")
		}
	df
}
