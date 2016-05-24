"gam.match" <-
function(x)
{
	if(is.list(x)) {
		junk <- Recall(x[[1]])
		if((nvar <- length(x)) == 1)
			return(list(o = junk$o, nef = junk$nef))
		else {
			o <- matrix(junk$o, length(junk$o), nvar)
			nef <- rep(junk$nef, nvar)
			for(i in 2:nvar) {
				junk <- Recall(x[[i]])
				o[, i] <- junk$o
				nef[i] <- junk$nef
			}
			names(nef) <- nn <- names(x)
			dimnames(o) <- list(NULL, nn)
			return(list(o = o, nef = nef))
		}
	}
	if(is.matrix(x)) {
		ats <- attributes(x)
		a <- ats$NAs
		ncols <- ats$ncols
		d <- dim(x)
		if(is.null(ncols))
			ncols <- d[2]
		if(ncols == 1)
			return(Recall(structure(x[, 1, drop = TRUE], NAs = a)))
		if(is.null(a)) {
			o <- seq(d[1])
			nef <- d[1]
		}
		else {
			nef <- d[1] - length(a)
			o <- rep(nef + 1, d[1])
			o[ - a] <- seq(nef)
		}
		return(list(o = as.integer(o), nef = as.integer(nef)))
	}
	else {
		a <- attributes(x)$NAs
		if(!is.null(a))
			x[a] <- NA
		xr <- signif(as.vector(x), 6)
		sx <- unique(sort(xr))
		nef <- as.integer(length(sx))
		if(nef <= 3)
			stop("A smoothing variable encountered with 3 or less unique values; at least 4 needed"
				)
		o <- match(xr, sx, nef + 1)
		o[is.na(o)] <- nef + 1
		return(list(o = as.integer(o), nef = as.integer(nef)))
	}
}
