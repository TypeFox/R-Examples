

flist <- function(f, feature) {
	method <- c("pstree", "cprob", "generate")
	stationary <- c(FALSE, TRUE, FALSE)
	cdata <- c(FALSE, FALSE, FALSE)
		
	fdb <- data.frame(method, stationary, cdata)

	return(fdb[fdb[,"method"]==f,feature])
}

