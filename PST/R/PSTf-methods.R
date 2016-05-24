## ======================
## number of observations
## ======================
setMethod("nobs", "PSTf", function(object) {
	
	res <- sum(object[[1]][["e"]]@n)
	
	return(res)
})

## summary
setMethod("summary", "PSTf", function(object, max.level) {

	if (missing(max.level)) { max.level <- length(object)-1 }

	stats <- PSTf.stats(object, max.level=max.level)
	res <- new("PST.summary",
		alphabet=object@alphabet,
		labels=object@labels,
		cpal=object@cpal,
		ns=as.integer(stats$ns),
		depth=as.integer(stats$depth),
		nodes=as.integer(stats$nodes),	
		leaves=as.integer(stats$leaves),
		freepar=as.integer((stats$nodes+stats$leaves)*(length(object@alphabet)-1))
	)

	return(res)
}
)

## Stats and summary
PSTf.stats <- function(PST, max.level) {
	stats <- list(ns=as.integer(0), leaves=as.integer(0), nodes=as.integer(0), depth=as.integer(0))

	stats$ns <- sum(PST[[1]][["e"]]@n)
	for (i in (max.level+1):1) {
		pruned.nodes <- pruned.nodes(PST[[i]])
		if (any(pruned.nodes)) {
			PST[[i]] <- PST[[i]][!pruned.nodes]
			new.leaves <- !names(PST[[i-1]]) %in% lapply(PST[[i]], node.parent) 
			PST[[i-1]][new.leaves] <- lapply(PST[[i-1]][new.leaves], node.leaf)
		}

		if (length(PST[[i]])>0) {
			if (i==(max.level+1)) {
				PST[[i]] <- lapply(PST[[i]], node.leaf)
			}

			if (stats$depth==0) { stats$depth <- i-1 }
			stats$nodes <- stats$nodes+sum(nodes.count(PST[[i]]))
			stats$leaves <- stats$leaves+sum(leaves.count(PST[[i]]))
		}
	}

	return(stats)
}


## ============================================================
## plot and print method builds a recursive version of x (PSTr)
## and call the corresponding method for class PSTr

## plot method
setMethod("plot", "PSTf", function (x, y=missing, max.level=NULL,
	nodePar = list(), edgePar = list(), 
	axis=FALSE, xlab = NA, ylab = if (axis) { "L" } else {NA}, 
	horiz = FALSE,  xlim=NULL, ylim=NULL, 
	withlegend=TRUE, ltext=NULL, cex.legend=1, use.layout=withlegend!=FALSE, legend.prop=NA, ...) {

	if (nrow(x@cdata)>0) {
		ccol <- cpal(x@cdata)
		cnames <- alphabet(x@cdata)

		if (attr(x@cdata,"nr") %in% names(x[[2]])) {
			ccol <- c(ccol, attr(x@cdata,"missing.color"))
			cnames <- c(cnames, attr(x@cdata,"nr"))
		}
		names(ccol) <- cnames

		if (!"stcol" %in% names(edgePar)) { edgePar[["stcol"]] <- ccol }
		if (!"c.cpal" %in% names(nodePar)) { nodePar[["c.cpal"]] <- ccol }
	}

	x <- as.pstree(x, max.level=max.level)

	plot(x, y=missing, max.level=max.level, nodePar=nodePar, edgePar=edgePar, axis=axis, 
		xlab=xlab, ylab=ylab, horiz=horiz, 
		xlim=xlim, ylim=ylim, 
		withlegend=withlegend, ltext=ltext, cex.legend=cex.legend, ...)

}
)

## print method
setMethod("print", "PSTf", function (x, max.level = NULL, ...) {

	x <- as.pstree(x, max.level=max.level)
	print(x, max.level = max.level, ...)
}
)

## node names
setMethod("nodenames", "PSTf", function (object, L) {
	
	if (missing(L)) {
		res <- NULL
		for (L in 1:length(object))  {
			res <- c(res, names(object[[L]]))
		}
	} else {
		res <- names(object[[L+1]])
	}
	
	return(res)
}
)




