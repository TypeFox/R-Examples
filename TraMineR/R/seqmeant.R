## ==============
## Mean durations
## ==============

seqmeant <- function(seqdata, weighted=TRUE, with.missing=FALSE, prop=FALSE, serr=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop("seqmeant: seqdata is not a sequence object, use seqdef function to create one")

	istatd <- suppressMessages(seqistatd(seqdata, with.missing=with.missing, prop=prop))

	weights <- attr(seqdata, "weights")

	if (!weighted || is.null(weights))
		weights <- rep(1, nrow(seqdata))
	## Also takes into account that in unweighted sequence objects created with
	## older TraMineR versions the weights attribute is a vector of 1
	## instead of NULL
	if (all(weights==1))
		weighted <- FALSE
	
	wtot <- sum(weights)

	mtime <- apply(istatd*weights,2,sum)
	
	res <- mtime/wtot

	res <- as.matrix(res)
	colnames(res) <- "Mean"

	col <- cpal(seqdata)
	if (with.missing) {
		col <- c(col, attr(seqdata,"missing.color"))
	}

	if(serr){
	  w2tot <- sum(weights^2)
	  vcent <- t(t(istatd) - mtime/wtot)
	  var <- apply(weights*(vcent^2),2,sum) * wtot/(wtot^2 - w2tot)
	  sd <- sqrt(var)
	  SE <- sqrt(var/wtot)
	  res <- cbind(res,var,sd,SE)
	  colnames(res) <- c("Mean", "Var", "Stdev", "SE")
	}

	attr(res,"nbseq") <- sum(weights)
	attr(res,"cpal") <- col
	attr(res,"xtlab") <- colnames(seqdata)
	attr(res,"weighted") <- weighted
	attr(res,"se") <- serr
	
	class(res) <- c("stslist.meant", "matrix")

	return(res)
}
