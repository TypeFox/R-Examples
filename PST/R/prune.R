## ==================================
## Pruning Probabilistic Suffix Trees
## ==================================

setMethod("prune", "PSTf", function(object, nmin, L, gain, C, keep, drop, state, delete=TRUE, lik=TRUE) {
	
	data <- object@data
	cdata <- object@cdata
	A <- alphabet(object)
	cpal <- cpal(object)
	labels <- stlab(object)
	segmented <- object@segmented
	group <- object@group

	if (!missing(gain) & missing(C)) {
		stop(" [!] please provide a cutoff value")
	} else if (missing(gain) & !missing(C)) {
		stop(" [!] missing 'gain' argument")
	}

	if (!missing(keep)) {
		if (has.cdata(object)) {
			c.A <- alphabet(object@cdata)
		} else {
			c.A <- object@alphabet
		}

		if (!inherits(keep,"stslist")) {
			keep <- seqdef(keep, alphabet=c.A, nr="#")
		}
		keep.sl <- seqlength(keep)
	}

	object <- as(object, "list")
	cnodes <- NULL

	message(" [>] pruning results: ")
	message("   ", format("[L]", width=5, justify="right"), format("[nodes]",width=9, justify="right"), 
		format("[pruned]", width=10, justify="right"))

	if (!missing(gain) && is.character(gain)) {
		if (gain=="G1") { gain <- G1 } else if (gain=="G2") { gain <- G2 }
	}

	for (i in length(object):2) {
		nodes <- object[[i]]
		parents <- object[[i-1]]
		## nodes <- lapply(nodes, function(x) {x@pruned[] <- FALSE; x})
		nbnodes <- unlist(lapply(nodes, function(x) { sum(!x@pruned) }))

		if (!missing(L) && i>(L+1) ) {
			nodes <- lapply(nodes, node.prune)
		} else {
			if (!missing(keep)) {
				if ( (i-1)>max(keep.sl) ) {
					nodes <- lapply(nodes, node.prune)
				} else {
					keep.tmp <- keep[keep.sl==i-1,, drop=FALSE]
					keep.list <- seqconc(keep.tmp)
					nodes <- lapply(nodes, node.keep, keep.list, clist=cnodes)
				}
			}
			if (!missing(state)) {
				state.tmp <- seqdecomp(names(nodes))
				state.tmp <- which(rowSums(state.tmp==state)>0)
				nodes[state.tmp] <- lapply(nodes[state.tmp], node.prune)
			}
			if (!missing(nmin)) {
				nodes <- lapply(nodes, node.nmin, nmin)
			}
			if (!missing(gain)) {
				nodes <- lapply(nodes, node.gain, plist=parents, gain=gain, C=C, clist=cnodes)
			}
		}

		pruned <- unlist(lapply(nodes, function(x) { sum(x@pruned) } ))
		plabel <- if (segmented) { " node segment(s) pruned" } else { " node(s) pruned" }

		message("   ", format(i-1, width=5), format(sum(nbnodes),width=9), format(sum(pruned), width=10))

		if (sum(pruned)>0) {
			cnodes <- lapply(nodes, delete.pruned)
			## Removing empty nodes
			pruned.id <- which(unlist(lapply(cnodes, function(x) { nrow(x@prob)==0 } )))
			if (length(pruned.id)>0) { cnodes <- cnodes[-pruned.id] }

			## 
			if (delete) {
				nodes <- cnodes
				cnodes <- NULL

				if (length(nodes)>0) {
					## Node segments having no more childrens set as leaves
					remaining <- lapply(nodes, function(x) { rownames(x@prob) })
					rplist <- unlist(lapply(nodes, node.parent))

					parents <- lapply(parents, set.leaves, remaining, rplist)
				} else {
					parents <- lapply(parents, function(x) { x@leaf[] <- TRUE; x })
				}

				## pnames <- unlist(lapply(nodes, node.parent))
				## parents <- lapply(parents, function(x, pnames) {if (!x@path %in% pnames) {x@leaf[] <- TRUE}; x}, pnames)
			} 
		} else if (!delete) {
				cnodes <- nodes
		}

		## =====================
		if (length(nodes)==0) {
			object <- object[-i]
		} else {
			object[[i]] <- nodes
		}
		object[[i-1]] <- parents
	}

	object <- new("PSTf", object, data=data, cdata=cdata, alphabet=A, cpal=cpal, labels=labels, 
		segmented=segmented, group=group, call=match.call(), logLik=as.numeric(NULL))

	## Loglik
	## Loglik
	if (lik) { object@logLik <- likelihood(object, log=TRUE) }

	return(object)
}
)

