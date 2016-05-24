## =============================================
## Computing Lth order conditional probabilities
## =============================================

setMethod("cprob", signature=c(object="stslist"), 
 	def=function(object, L, cdata=NULL, context, stationary=TRUE, nmin=1, prob=TRUE, weighted=TRUE, 
		with.missing=FALSE, to.list=FALSE) {

	debut <- Sys.time()

	if (!is.null(cdata) & !flist("cprob", "cdata")) {
			stop(" [!] argument cdata not available", call.=FALSE)
	}

	statl <- alphabet(object)
	## We see if the name of one of the states conflicts with names used to label the results 
	if (any(statl=="[n]")) { stop(" [!] found the symbol '[n]' in the alphabet, which conflicts with internal names used in PST. Please rename this state.") } 
	nr <- attr(object,"nr")
	if (with.missing) { statl <- c(statl, nr) }
	nbetat <- length(statl)
	sl <- seqlength(object)
	sl.max <- max(sl)
	if (L>(sl.max-1)) { stop(" [!] sequence length <= L")}
	nbseq <- nrow(object)

	## Weights
	weights <- attr(object, "weights")

	if (!weighted || is.null(weights)) {
		weights <- rep(1, nrow(object))
	}

	## Turning object into a matrix
	object <- as.matrix(object)

	if (!missing(context)) {
		tmp <- seqdecomp(context)
		if (any(!tmp %in% statl) & context!="e") {
			stop(" [!] one or more symbol in context not in alphabet")
		}
		L <- ncol(tmp)
	} 

	message(" [>] ", nbseq, " sequences, min/max length: ", min(sl), "/", max(sl))

	states <- factor(object[,(L+1):sl.max], levels=statl)
	contexts <- matrix(nrow=nbseq, ncol=sl.max-L)
	if (L==0) {
		contexts[] <- "e"
	} else {
		if (is.null(cdata)) { cdata <- object } else { cdata <- as.matrix(cdata) }

		for (p in (L+1):sl.max) {
			contexts[, p-L] <- cdata[, (p-L)]
			if (L>1) {
				for (c in (L-1):1) {
					contexts[, p-L] <- paste(contexts[, p-L], cdata[, (p-c)], sep="-")
				}
			}	
		}
		## This version using apply is slower 
		## for (p in (L+1):sl.max) {
		## 	contexts[, p-L] <-	apply(cdata[, (p-L):(p-1), drop=FALSE], 1, paste, collapse="-")
		## }
	}
	contexts <- as.vector(contexts)

	## inflating weight vector to match number of contexts
	weights <- rep(weights, ncol(object)-L)

	if (!missing(context)) {
		sel <- contexts==context
		contexts <- contexts[sel]
		states <- states[sel]
		weights <- weights[sel]
	} 

	message(" [>] computing prob., L=", L, ", ", length(unique(contexts)), " distinct context(s)") 

	if (stationary) {
			freq <- xtabs(weights ~ contexts+states)[,,drop=FALSE]
			if (prob) { freq <- freq/rowSums(freq) }
			n <- rowSums(xtabs( ~ contexts+states)[,,drop=FALSE])	
			res <- cbind(freq, n)
			colnames(res) <- c(statl, "[n]")

			if (nmin>1) {
				nmin.del <- which(res[,"[n]"]<nmin)
				if (length(nmin.del)>0) {
					res <- res[-nmin.del,,drop=FALSE]
					message(" [>] removing ", length(nmin.del), " context(s) where n<", nmin) 
				}
			}
	} else {
		t <- (L+1):sl.max
		pos <- matrix(t, ncol=length(t), nrow=nbseq, byrow=T)
		pos <- as.vector(pos)
		pos <- factor(pos)
		if (!missing(context)) { pos <- pos[sel] }

		tmat <- xtabs(weights ~ pos+states+contexts)
		n <- xtabs( ~ pos+states+contexts)
		context.list <- dimnames(tmat)$contexts

		## Creating a list with conditional prob+n for each context
		res <- lapply(context.list, 
			function(idx) {
				## If only one row tmat[,,idx] cannot be extracted as a matrix !!
				if (nrow(tmat[,,idx, drop=FALSE])==1) {
					freq <- t(as.matrix(tmat[,,idx]))
					rownames(freq) <- rownames(tmat[,,idx, drop=FALSE])
					fs <- sum(n[,,idx])
				} else {
					freq <- tmat[,,idx]
					fs <- rowSums(n[,,idx])
				}

				if (prob) { freq <- freq/rowSums(freq) }

				pplist <- cbind(freq, "[n]"=fs)
				## colnames(pplist) <- c(statl, "n")
				## Sorting by position
				pplist <- pplist[order(as.numeric(rownames(pplist))),, drop=FALSE]

				nmin.del <- which(pplist[,"[n]"]<nmin)
				if (length(nmin.del)>0) {
					pplist <- pplist[-nmin.del,,drop=FALSE]
				}

				return(pplist)
			}
		)

		names(res) <- context.list
	
		if (L>0) {
			## eliminating empty elements
			empty <- which(unlist(lapply(res, function(x) { is.null(x) || nrow(x)==0 } )))
			if (length(empty)>0) { res <- res[-empty] }
		}
	}

	## eliminating patterns containing missing states if with.missing=FALSE
	if (L>0 & !with.missing) { 
			## if nr is a 'grep' special character 
			if (nr %in% c("?", "*")) { nr <- paste("\\",nr, sep="") }
			hasMiss <- if (stationary) { grep(nr, rownames(res)) } else { grep(nr, names(res)) }
			if (length(hasMiss)>0) { 
				message(" [>] removing ", length(hasMiss), " context(s) containing missing values") 
				res <- if (stationary) { res[-hasMiss, ,drop=FALSE] } else { res[-hasMiss] }
			}
	}

	fin <- Sys.time()
	message(" [>] total time: ", format(round(fin-debut, 3)))

	if (stationary & to.list) {
		res <- lapply(1:nrow(res), function(i) res[i,,drop=FALSE])
		nodes.names <- unlist(lapply(res, rownames))
		names(res) <- nodes.names
		res <- lapply(res, function(x) {rownames(x) <- NA; x} )
	}

	return(res)
}
)
