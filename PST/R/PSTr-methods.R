## Initialize method

setMethod("initialize", "PSTr", function(.Object, path, counts, n, order, ymin, index, ...) {

	idx <- NULL

	if (all(is.na(index[,"group"])) & all(is.na(index[,"position"]))) { 
		idx <- "S1" 
	} else {
		if (!all(is.na(index[,"group"]))) { idx <- paste("G", index[,"group"], sep="") }
		if (!all(is.na(index[,"position"]))) { idx <- paste(idx, paste("P", index[,"position"], sep=""), sep="") }
	} 
	
	rownames(index) <- idx
	.Object@index <- index
	
	nd <- nrow(counts)
	.Object@leaf <- matrix(FALSE, nrow=nd, dimnames=list(idx, ""))
	.Object@pruned <- matrix(FALSE, nrow=nd, dimnames=list(idx, ""))
	rownames(counts) <- idx
	rownames(n) <- idx

	counts <- as.matrix(counts)
	.Object@prob <- counts/rowSums(counts)

	## Smoothing
	if (!is.null(ymin) & any(.Object@prob==0)) {
		ids <- which(rowSums(.Object@prob==0)>0)
		.Object@prob[ids,] <- ((1-(ncol(counts)*ymin)) * .Object@prob[ids,]) + ymin
	}
	callNextMethod(.Object, path=path, counts=counts, n=n, order=order, ...)
}
)

## TAKEN FROM dengrogram method
## The ``generic'' method for "[["  (identical to e.g., "[[.POSIXct"):
## --> subbranches are pstrees as well!
setMethod("[[", "PSTr", function(x, i, drop = TRUE, root.attr=FALSE) {
	## cl <- class(x)
	alphabet <- x@alphabet
	cpal <- x@cpal
	labels <- x@labels

	x <- unclass(x)
	x <- x[[i]]

	if (!is.null(x) & root.attr) {
		x@alphabet <- alphabet
		x@cpal <- cpal
		x@labels <- labels
	}

	x
}
)

## Summary
setMethod("summary", "PSTr", function(object, max.level=NULL, segmented=TRUE) {

	stats <- PSTr.stats(object, max.level=max.level, segmented=segmented)

	res <- new("PST.summary",
		alphabet=object@alphabet,
		labels=object@labels,
		cpal=object@cpal,
		ns=as.integer(object@n),
		depth=as.integer(stats$depth),
		nodes=as.integer(stats$nodes),	
		leaves=as.integer(stats$leaves),
		freepar=as.integer((stats$nodes+stats$leaves)*(length(object@alphabet)-1))
	)

	return(res)
}
)

PSTr.stats <- function(PST, max.level, segmented=TRUE) {
	stats <- list(leaves=as.integer(0), nodes=as.integer(0), depth=as.integer(0))

	childrens <- which.child(PST)

	stats$depth <- attr(PST,"order")

	if (!is.null(max.level) && PST@order==max.level) {
		PST@leaf[] <- TRUE
	}

	if (segmented) {
		stats$leaves <- stats$leaves+sum(PST@leaf)
	} else {
		stats$leaves <- stats$leaves+all(PST@leaf)
	}

	if (any(!PST@leaf)) {
		stats$depth <- stats$depth+1
		if (segmented) {
			stats$nodes <- stats$nodes+sum(!PST@leaf)
		} else {
			stats$nodes <- stats$nodes+1
		}

		for (i in childrens) {
			tmp <- PSTr.stats(PST[[i]], max.level=max.level, segmented=segmented)
			stats$leaves <- stats$leaves+tmp$leaves
			stats$nodes <- stats$nodes+tmp$nodes
			if (tmp$depth>stats$depth) (stats$depth  <- tmp$depth)
		}
	}
	
	return(stats)
}


setMethod("show", "PST.summary", function(object) {
	alphabet <- object@alphabet
	nbstates <- length(alphabet)
	cpal <- object@cpal
	labels <- object@labels
	ns <- object@ns

	cat(" [>] alphabet (state labels): ","\n")
	maxstatedisplay <- 12
	for (i in 1:min(nbstates,maxstatedisplay))
		cat("     ",i, "=", object@alphabet[i], " (", labels[i], ")","\n", sep="")
	if (nbstates>12) message("      ...")
	cat(" [>] PST fitted on", ns, "symbols \n") 
	cat(" [>] max. depth:", object@depth,"\n")
	cat(" [>] ", object@nodes, " internal node(s), ", object@leaves, " leave(s)\n", sep="")
	cat(" [>] (", object@nodes,"+",object@leaves,") * (",nbstates, "-1) = ", object@freepar, " free parameters \n", sep="")
}
)

which.child <- function(PST) {
	class <- class(PST)
	child.list <- names(PST)[unlist(lapply(PST, is, class))]
	return(child.list)
}

