## BUILDING A PROBABILISTIC SUFFIX TREE

setMethod("pstree", signature="stslist", 
	function(object, group, L, cdata=NULL, stationary=TRUE, nmin=1, ymin=NULL, weighted=TRUE, with.missing=FALSE, lik=TRUE) {

	debut <- Sys.time()

	if (!stationary & !flist("pstree", "stationary")) {
			stop(" [!] 'stationary=FALSE' not implemented", call.=FALSE)
	}

	if (!is.null(cdata) & !flist("pstree", "cdata")) {
			stop(" [!] argument cdata not allowed", call.=FALSE)
	}

	if (missing(L)) { L <- max(seqlength(object))-1 }

	if (!weighted || is.null(attr(object, "weights"))) { attr(object, "weights") <- rep(1.0, nrow(object)) }
	## Also takes into account that in unweighted sequence objects created with 
	## older TraMineR versions the weights attribute is a vector of 1
	## instead of NULL  
	if (all(attr(object, "weights")==1)) { weighted <- FALSE }

	A <- alphabet(object)
	if (any(A=="[n]")) { stop(" [!] found the symbol '[n]' in the alphabet, which conflicts with internal names used in PST. Please rename this state.") } 
	StCol <- cpal(object)
	StLab <- stlab(object)
	sl <- seqlength(object)

	if (with.missing) { 
		A <- c(A, attr(object, "nr"))
		StCol <- c(StCol, attr(object, "missing.color"))
		StLab <- c(StLab, "missing")
	}
	names(StCol) <- A

	## Segmented PST
	if (!missing(group)) {
		if (!is.factor(group)) { group <- factor(group) }
		nbgroup <- length(levels(group))
		segmented <- TRUE
		message(" [>] ", nbgroup, " groups")
	} else {
		group <- factor(NULL)
		segmented <- FALSE
	}

	## 
	nodes.list <- vector("list", length=L+1)
	
	message(" [>] ", nrow(object), " sequence(s) - min/max length: ", min(sl),"/",max(sl))
	message(" [>] max. depth L=", L, ", nmin=", nmin, if (!is.null(ymin)) { paste(", ymin=", ymin, sep="") })
	message("   ", format("[L]", width=5, justify="right"), format("[nodes]",width=9, justify="right"))

	for (i in 0:L) {
		if (segmented) {
			tmp <- NULL
			idx <- NULL
			for (g in 1:nbgroup) {
				data <- object[group==levels(group)[g],]
				tmp.cdata <- if (!is.null(cdata)) { cdata[group==levels(group)[1],] } else { NULL }
				ccounts <- suppressMessages(cprob(data, L=i, cdata=tmp.cdata, stationary=stationary,
					nmin=nmin, prob=FALSE, weighted=weighted, with.missing=with.missing, to.list=TRUE))

				ccounts <- lapply(ccounts, function(x) cbind(x, group=g, position=as.integer(rownames(x))))
				tmp <- merge.cprob(tmp, ccounts)
			}
		} else {
			tmp <- suppressMessages(cprob(object, L=i, cdata=cdata, stationary=stationary,
				nmin=nmin, prob=FALSE, weighted=weighted, with.missing=with.missing, to.list=TRUE))
			tmp <- lapply(tmp, function(x) { cbind(x, group=NA, position=as.integer(rownames(x))) } )
		}

		## 
		nodes.names <- names(tmp)
		tmp.list <- lapply(seq_len(length(tmp)), function(n) 
			new("PSTr", path=nodes.names[n], counts=tmp[[n]][,A, drop=FALSE], n=tmp[[n]][,"[n]", drop=FALSE], 
				order=i, ymin=ymin, index=tmp[[n]][,c("group","position"),drop=FALSE])
		)

		nbnodes <- length(tmp.list)
		names(tmp.list) <- nodes.names

		## Listing children path+position
		if (i>0) { 
			parents <- nodes.list[[i]]
			child.list <- lapply(tmp.list, function(x) { rownames(x@prob) } )
			rplist <- unlist(lapply(tmp.list, node.parent))
			
			## Nodes having no children set as leaves
			nodes.list[[i]] <- lapply(parents, set.leaves, child.list, rplist)
		}
			
		message("   ", format(i, width=5), format(nbnodes, width=9))
		nodes.list[[i+1]] <- tmp.list
	}
	
	## At max depth all nodes are leaves
	nodes.list[[L+1]] <- lapply(nodes.list[[L+1]], node.leaf)

	if (is.null(cdata)) { 
		cdata <- object[-(1:nrow(object)),] 
	} else {
		if (with.missing) {
			cdata <- seqdef(cdata, nr="#", alphabet=c(alphabet(cdata), attr(cdata, "nr")),
				labels=c(stlab(cdata), "missing"), cpal=c(cpal(cdata), attr(cdata, "missing.color")),
				xtstep=attr(cdata, "xtstep"))
		}
	} 

	res <- new("PSTf", nodes.list, data=object, cdata=cdata, alphabet=A, cpal=StCol, labels=StLab, 
		segmented=segmented, group=group, call=match.call(), logLik=as.numeric(NULL))

	## likelihood
	if (lik) { res@logLik <- likelihood(res, log=TRUE) }

	fin <- Sys.time()
	message(" [>] total time: ", format(round(fin-debut, 3)))

	return(res)
}
)

## 
merge.cprob <- function(x,y) {
		if (is.null(x)) {
			res <- y
		} else {
			for (i in names(y)) {
				if (i %in% names(x)) {
					x[[i]] <- rbind(x[[i]], y[[i]])
				} else {
					x[[i]] <- y[[i]]
				}
			}
			res <- x
		}
}






