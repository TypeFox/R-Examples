# ANOVA comparison of several fSRM objects; 
# Takes a list of lavaan models (can include NULLs), and returns the usual anova object
anovaList <- function(modellist) {
	# a0: list without NULLs
	
	mlist2 <- sapply(modellist, function(x) return(x$fit))
	
	mods <- mlist2[!sapply(mlist2, function(x) is.null(x))]
	mods <- mods[!sapply(mods, function(x) !inspect(x, "converged"))]
	
	if (length(mods) == 0) {
		return(list(n.mods=0))
	}
		
    # put them in order (using number of free parameters)
    nfreepar <- sapply(mods, function(x) x@Fit@npar)
    if(any(duplicated(nfreepar))) { ## FIXME: what to do here?
        # what, same number of free parameters?
        # maybe, we need to count number of constraints
        ncon <- sapply(mods, function(x) { nrow(x@Model@con.jac) })
        nfreepar <- nfreepar - ncon
    }

    mods <- mods[order(nfreepar, decreasing = TRUE)]		
		
	pStr <- sapply(1:length(mods), function(x){ 
		if(x==1) {
			paste("mods[[",x,"]]",sep = "")
		} else {
			paste("force(mods[[",x,"]])",sep = "")
		}
	})
	pStr2 <- paste0("anova(", paste(pStr, collapse=", "), ")")
	
	a1 <- eval(parse(text = pStr2))
	
	if (length(mods) > 1) {
		rownames(a1) <- names(modellist)
	}
	
	attr(a1, "n.mods") <- length(mods)
	return(a1)
}