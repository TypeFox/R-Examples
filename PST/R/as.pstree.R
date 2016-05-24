## ===========================================
## from "flat" to recursive PST representation
## ===========================================

as.pstree <- function(object, max.level=NULL) {

	if (is.null(max.level)) { max.level <- length(object)-1 }

	debut <- Sys.time()

	if ((length(object)-1)<max.level) { 
		max.level <- length(object)-1
	}

	message(" [>] building 'PSTr' representation, max. depth=", max.level, "...", appendLF=FALSE)

	N0 <- object[[1]][["e"]]
	A <- alphabet(object)
	N0@alphabet <- A
	N0@labels <- stlab(object)
	N0@cpal <- cpal(object)

	for (i in 2:(max.level+1)) {
		id.comp <- names(object[[i]])
		id <- seqdecomp(id.comp)
		nbseg <- length(id.comp)

		if (nbseg>0) {
			A.sort <- match(id[,1], A)
			node <- id[,ncol(id):1, drop=FALSE]
			node <- node[order(A.sort), ,drop=FALSE]
			id.comp <- id.comp[order(A.sort)]
				
			for (j in 1:nbseg) {
				if (getOption("verbose")==TRUE) { message(" [>] adding node: ", rev(node[j,]))}
					## The node's parent is the longest suffix of the string
				N0[[node[j,]]] <- object[[i]][[id.comp[j]]]
			}
		}

	}

	fin <- Sys.time()
	message(" (", format(round(fin-debut, 3)), ")")

	return(N0)
}

##
setAs(from="PSTf", to="PSTr", def=function(from) as.pstree(from))


