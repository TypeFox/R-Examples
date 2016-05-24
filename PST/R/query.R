## Extracting the probability of observing each symbol in the alphabet after a given subsequence (context)

setMethod("query", signature=c(object="PSTf"), 
	def=function(object, context, state, output="prob", exact=FALSE) {
		A <- attr(object, "alphabet")
		cA <- if (nrow(object@cdata)>0) { alphabet(object@cdata) } else { A }

		if (missing(context) || context=="e") {
			context <- "e"
			idxl <- 1
		} else {
			sd <- unlist(strsplit(context, split="-"))
      			context <- paste(sd, collapse="-")
			idxl <- length(sd)+1
			if (any(!sd %in% cA)) {
				stop(" [!] one or more symbol not in alphabet")
			} 
		}

		if (exact && (!context %in% names(object[[idxl]]) || object[[idxl]][[context]]@pruned)) {
			message( "[>] node is not in the tree")
			res <- NULL
		} else {
			if (idxl>length(object)) { 
				idxl <- length(object) 
				sd <- sd[(length(sd)-(idxl-2)):length(sd)]
        			context <- paste(sd, collapse="-")
			}

			while (!context %in% names(object[[idxl]]) || object[[idxl]][[context]]@pruned) {
				idxl <- idxl-1
        			sd <- sd[2:length(sd)]
				context <- if (idxl>1) { paste(sd, collapse="-") } else {"e"}
			}

			res <- object[[idxl]][[context]]

			if (output=="prob") {
				res <- res@prob
			} else if (output=="counts") {
				res <- res@counts
			} else if (output=="n") {
				res <- res@n
			}
			message(" [>] retrieving from node: ", paste(context, collapse="-"))
		}

		if (output %in% c("prob", "counts") & !is.null(res)) {
			if (!missing(state)) {
				res <- res[, which(A==state), drop=FALSE]
			}
			res <- new("cprobd", res, context=context)
		}

		return(res)
	}
)




