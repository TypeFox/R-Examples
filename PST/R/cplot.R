## Probability distribution of a node

setMethod("cplot", signature="PSTf", 
	def=function(object, context, state, main=NULL, all=FALSE, x.by=1, y.by=0.2, by.state=FALSE, ...) {

		A <- object@alphabet
		cpal <- c(object@cpal)
		oolist <- list(...)
		
		if (!missing(state) && is.character(state)) { 
			state <- which(state==A)
		}

		node <- query(object, context, output="all")
		prob <- node@prob
		pruned <- node@pruned

		if (all) {
			seglist <- object[[1]][[1]]@index
		} else {
			seglist <- node@index
		}

		if (is.null(main)) {
			main <- paste("Node ", node@path)
		}

		plot(NULL, axes=FALSE, ylab="Prob", xlim=c(1,(nrow(seglist)+1)), ylim=c(1,0), main=main, ...)

		## prob matrix is reversed because we are using the function for plotting the tree nodes (yaxis is reversed) 
		plotNodeProb(1, 1, nrow(seglist)+1, 0, prob=prob, seglist=seglist, cpal=cpal, pruned=pruned, 
			by.state=by.state, index=node@index)

		if (length(seglist)>1) { 
			slab.pos <- seq(1, nrow(seglist), by=x.by)
			axis(1, at=slab.pos+0.5, labels=seglist[slab.pos]) 
		}
		axis(2, at=seq(0,1, by=y.by), labels=rev(seq(0,1, by=y.by)))
	}
)

