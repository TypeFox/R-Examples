search.match <- function(data, treatment, vars, depth=3, min.vars =1, group=1, useCP, ...)
{
	#if(!require(MatchIt))
	# stop("You need MatchIt package to use this tool")
	
	if(missing(data) | missing(treatment) | missing(useCP))
		stop("One or more of these arguments are missing: treatment, data, useCP.")

	xvars <- NULL

	if(missing(vars)){
		idx <- which(vars == treatment)
		if(length(idx)==0)
			stop("Treatment variable wrong o missing")
		xvars <- colnames(data)[-idx]
	} else {
		xvars <- unique(vars)
	}  

	start <- paste(treatment,"~") #conditional effects

	total <- 1
	for (i in min.vars:(length(xvars)-1))
		total <- total + ncol(combn(xvars, i)) 
 
	for (i in 2:depth)
		total <- total + ncol(combn(xvars, i)) 
	
	
	c.idx <- NULL	
	for(i in xvars){
		if(is.numeric(data[[i]]) | is.integer(data[[i]]))
			c.idx <- c(c.idx, i)
	}

	n.cvars <- length(c.idx)
	cvars <- NULL
	if(n.cvars>0){	
		cvars <- c.idx
		for (i in 1:min(depth, n.cvars))
			total <- total + ncol(combn(cvars, i)) 
	}
	
	cat("\n I'm going to run", total, "different matching solutions!\n\n")
	pb <- txtProgressBar(min = 1, max = total, initial = 1, style = 3)

	L1 <- rep(as.numeric(NA), total)
	n <- rep(as.numeric(NA), total)
	frml <- character(total)
	counter <- 1
	for (i in min.vars:(length(xvars)-1)) {
		allsubset <- combn(xvars, i)
		for (j in 1:ncol(allsubset)) {
			ftmp <- start
			for (k in 1:i)
			ftmp <- paste(ftmp, "+", allsubset[k,j])
			frml[counter] <- ftmp
			ftmp <- as.formula(ftmp)
			setTxtProgressBar(pb, counter)
			psm <- try(matchit(ftmp, data=data, ...), silent = TRUE)
			if (class(psm) != "try-error") {
				L1[counter] <- L1.meas(data[[treatment]], data[xvars], breaks=useCP, weights=psm$weights)$L1
				n[counter] <- sum(psm$weights>0 & data[[treatment]]==group)
			}
			counter <- counter + 1
		}
	}
	
	
	mfull <- paste("treated","~",paste(xvars,collapse=" + "))

	for (i in 2:3) {
		allsubset <- combn(xvars, i)
		for (j in 1:ncol(allsubset)) {
			ftmp <- paste( mfull , "+", paste(allsubset[,j],collapse="*"))
			frml[counter] <- ftmp
			ftmp <- as.formula(ftmp)
			setTxtProgressBar(pb, counter)
			psm <- try(matchit(ftmp, data=data, ...), silent = TRUE)
			if (class(psm) != "try-error") {
				L1[counter] <- L1.meas(data[[treatment]], data[xvars], breaks=useCP, weights=psm$weights)$L1
				n[counter] <- sum(psm$weights>0 & data[[treatment]]==group)
			}
			counter <- counter + 1
		}
	}
	
	if(n.cvars>1){
		for (i in 2:min(3, n.cvars)) {
			allsubset <- combn(cvars, i)
			for (j in 1:ncol(allsubset)) {
				ftmp <- paste( mfull , "+", paste(sprintf("+ %s^%d",allsubset[,j],2), collapse=""))
				frml[counter] <- ftmp
				ftmp <- as.formula(ftmp)
				setTxtProgressBar(pb, counter)
				psm <- try(matchit(ftmp, data=data, ...), silent = TRUE)
				if (class(psm) != "try-error") {
					L1[counter] <- L1.meas(data[[treatment]], data[xvars], breaks=useCP, weights=psm$weights)$L1
					n[counter] <- sum(psm$weights>0 & data[[treatment]]==group)
				}
				counter <- counter + 1
			}
		}
	}
	
	
	close(pb)

	return( list(n=n, L1=L1, frml=frml, CP=useCP) ) 
}



