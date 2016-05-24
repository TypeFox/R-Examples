## Computing sequence probability

setMethod("predict", signature=c(object="PSTf"), 
	def=function(object, data, cdata, group, L=NULL, p1=NULL, output="prob", decomp=FALSE, base=2) {

	mlist <- c("prob", "SIMn", "logloss", "SIMo")

	if (!output %in% mlist) {
		stop(" output must be one of: ", mlist)
	}

	if (object@segmented) {
		sl <- seqlength(data)
		nbgroup <- length(levels(object@group))
		
		if (missing(group)) {
			prob <- matrix(nrow=nrow(data), ncol=nbgroup)
			colnames(prob) <- levels(object@group)
			rownames(prob) <- rownames(data)
			message(" [>] cross prediction for ", nrow(data), " sequence(s) - ", nbgroup, " models")
		} else {
			group <- factor(group)
			if (length(group)!=nrow(data)) {
				stop(" group must contain one value for each sequence in data")
			}
			if (decomp) {
				prob <- matrix(nrow=nrow(data), ncol=max(sl))
				colnames(prob) <- colnames(data)
			} else {
				prob <- matrix(nrow=nrow(data), ncol=1)
				colnames(prob) <- output
			}
			rownames(prob) <- rownames(data)
		}

		for (g in 1:nbgroup) {
			pst <- subtree(object, group=g)

			if (!missing(group)) {
				message(" [>] group ", g, ": ",appendLF=FALSE)
				group.idx <- which(group==levels(group)[g])
				prob[group.idx,] <- predict(pst, data[group.idx,], L=L, p1=p1, output=output, 
					decomp=decomp, base=base)
			} else {
				message(" [>] predicting with submodel ", g, ": ",appendLF=TRUE)
				prob[,g] <- predict(pst, data, L=L, p1=p1, output=output, decomp=decomp, base=base)
			}
		}
	
		return(prob)
	} else {
		A <- object@alphabet
		n <- nrow(data)
		sl <- seqlength(data)

		message(" [>] ", n, " sequence(s) - min/max length: ", min(sl),"/",max(sl))

		if (any(as.data.frame(data)==attr(data,"nr"))) {
			message(" [>] found missing value in sequence data")
			if (!attr(data, "nr") %in% A) {
				message(" [>] PST built without missing value as additional state")
				if (!decomp) {
					message("    [>] sequence probabilities not calculated on the same number of values")
				}
			}
		}

		debut <- Sys.time()


		if (min(sl)!=max(sl) & output=="prob" & !decomp) {
			message(" [!] sequences have unequal lengths")
		}

		stationary <- if (is.stationary(object)) { TRUE } else { FALSE }

		if (is.null(L)) { L <- length(object)-1 } 
		message( " [>] max. context length: L=", L)

		if (output=="SIMo") { P0 <- predict(object, data, L=0, output="prob", decomp=TRUE) }

		## getting contexts of max length L for each state
		contexts <- if (!missing(cdata) && nrow(cdata)>0) { context(object@cdata, L=L) } else { context(data, L=L) }
		states <- as.matrix(data)

		if (stationary) {
			context.table <- unlist(lapply(object[1:(L+1)], names))
			pruned.nodes <- unlist(lapply(object[1:(L+1)], pruned.nodes))
			if (any(pruned.nodes)) { 
				message(" [>] ", sum(pruned.nodes), "pruned nodes removed")
				context.table <- context.table[!pruned.nodes] 
			}

			## getting contexts of max length L for each state
			contexts <- as.vector(contexts)
			states <- as.vector(states)
			prob <- vector("numeric", length=length(states))
			prob[] <- NA
			unique.contexts <- unique(contexts)

			context.idx <- match(contexts, unique.contexts)

			## taking longest suffix until context found in PST
			unmatched <- !unique.contexts %in% context.table
			while (sum(unmatched>0)) {
				tmp <- seqdecomp(unique.contexts[unmatched])
				if (ncol(tmp)>1) {
					unique.contexts[unmatched] <- seqconc(tmp[,2:ncol(tmp), drop=FALSE])
					unique.contexts[unique.contexts==""] <- "e"
				} else { 
					unique.contexts[unmatched] <- "e"
				}
	
				## we may have reduced the number of distinct contexts
				tmp <- unique(unique.contexts)
				unique.match <- match(unique.contexts, tmp)
				context.idx <- unique.match[context.idx]

				##
				unique.contexts <- tmp
				unmatched <- !unique.contexts %in% context.table
				## print(unique.contexts[unmatched])
			}
	
			message(" [>] found ", length(unique.contexts), " distinct context(s)")
  
			for (p in 1:length(unique.contexts)) {
				context <- unique.contexts[p]
				context.eq <- which(context.idx==p)
      
				if (context=="e") {
					if (!is.null(p1)) {
						tmp <- p1
					} else {
						tmp <- object[[1]][["e"]]@prob[1,]
					}
				} else {
					sd <- unlist(strsplit(context, split="-"))
					idxl <- length(sd)+1	

					tmp <- object[[idxl]][[context]]@prob
				}

				tmp <- as.numeric(tmp)

				for (s in 1:length(tmp)) {
        			state.eq <- states[context.eq]==A[s]
					prob[context.eq][state.eq] <- tmp[s]
				}
			}
  			prob <- matrix(prob, ncol=max(sl))
		} else {
			prob <- matrix(nrow=nrow(data), ncol=max(sl))

			for (p in 1:max(sl)) {
				context.table <- NULL

				message(" [>] position ", p, " ...", appendLF=FALSE)

				for (i in 1:min(p,(L+1))) {
					tmp.plist <- unlist(lapply(object[[i]], function(x) p %in% x@index[,"position"]))
					tmp.pplist <- names(tmp.plist)[tmp.plist]
					context.table <- c(context.table, tmp.pplist)
				}

				unique.contexts <- unique(contexts[,p])

				contexts.idx <- match(contexts[,p], unique.contexts)

				## taking longest suffix until prefix found in PST
				unmatched <- !unique.contexts %in% context.table
				while (sum(unmatched>0)) {
					tmp <- seqdecomp(unique.contexts[unmatched])
					if (ncol(tmp)>1) {
						unique.contexts[unmatched] <- seqconc(tmp[,2:ncol(tmp), drop=FALSE])
						unique.contexts[unique.contexts==""] <- "e"
					} else { 
						unique.contexts[unmatched] <- "e"
					}

					## we may have reduced the number of distinct contexts
					tmp <- unique(unique.contexts)
					unique.match <- match(unique.contexts, tmp)
					contexts.idx <- unique.match[contexts.idx]

					##
					unique.contexts <- tmp
					unmatched <- !unique.contexts %in% context.table
				}
	
				message(" extracting ", length(unique.contexts), " context(s)")

				for (j in 1:length(unique.contexts)) {
					prefix <- unique.contexts[j]
					prefix.eq <- which(contexts.idx==j)
      
					if (prefix=="e") {
						if (!is.null(p1)) {
							tmp <- p1
						} else {
							tmp <- object[[1]][["e"]]@prob[p,]
						}
					} else {
						sd <- unlist(strsplit(prefix, split="-"))
						idxl <- length(sd)+1
						node <- object[[idxl]][[prefix]]
						idxp <- which(node@index[,"position"]==p)

						tmp <- node@prob[idxp,]
					}

					tmp <- as.numeric(tmp)

					for (s in 1:length(tmp)) {
        				state.eq <- states[prefix.eq, p]==A[s]
						prob[prefix.eq, p][state.eq] <- tmp[s]
					}
				}
			}
		}

		## Preparing final matrix
		rownames(prob) <- rownames(data)
		colnames(prob) <- colnames(data)

		if (output=="logloss") {
			prob <- -log(prob, base=base)
		} else if (output=="SIMn") {
			prob <- log(prob, base=base)
		} else if (output=="SIMo") {
			prob <- log(prob/P0, base=base)
		}

		if (!decomp) {
			if (output=="prob") {
				prob <- apply(prob,1, rowProds)
			} else if (output %in% c("logloss", "SIMn")) {
				prob <- rowSums(prob, na.rm=TRUE)/rowSums(!is.na(prob))
			} else if (output=="SIMo") {
				prob <- rowSums(prob, na.rm=TRUE)
			}

			## if only one sequences we return a matrix as well 
			if (is.null(dim(prob))) { prob <- matrix(prob, nrow=nrow(data)) }
			rownames(prob) <- rownames(data)
			colnames(prob) <- output
		}

		fin <- Sys.time()
		message(" [>] total time: ", format(round(fin-debut, 3)))

		return(prob)
	}
}
)


rowProds <- function(x) {

	vpos <- which(!is.na(x))

	if (length(vpos)>0) {
		p <- 1		
		for (i in vpos) {
			p <- p*x[i]
		}
	} else {
		p <- NULL
	}

	return(p)
}








