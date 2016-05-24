## Model selection according to the pruning cutoff 

setMethod("tune", signature=c(object="PSTf"), 
	def=function(object, gain="G2", C, criterion="AIC", output="PST") {

	if (!inherits(object, "PSTf") || missing(object)) {
		stop(" [!] please provide the name of a valid PST object", call.=FALSE)
	}

	debut <- Sys.time()
	nbmod <- length(C)

	nodes <- vector(mode="integer", length=nbmod) 
	leaves <- vector(mode="integer", length=nbmod) 
	k <- vector(mode="integer", length=nbmod) 
	IC <- vector(mode="numeric", length=nbmod)
	IC[] <- NA

	C <- sort(C)
	IC.function <- if (criterion=="AICc") { "AIC" } else { criterion }

	for (i in 1:nbmod) {
		suppressMessages(pst <- prune(object, gain=gain, C=C[i]))
		pst.sum <- summary(pst)
		nodes[i] <- pst.sum@nodes
		leaves[i] <- pst.sum@leaves
		k[i] <- pst.sum@freepar
		n <- pst.sum@ns

		IC.args <- list(object=pst)

		suppressMessages(pst.IC <- do.call(IC.function, args=IC.args))

		if (criterion=="AICc") {
			pst.IC <- pst.IC+((2*k[i]*(k[i]+1))/(n-k[i]-1))
		}

		IC[i] <- pst.IC
		
		message(" [>] model ",i, ": ", criterion,"=", round(pst.IC,2), " (C=", round(C[i],2),")")

		## 
		if (pst.IC==min(IC, na.rm=TRUE)) {
			id.best <- i
			pst.best <- pst
		}
	}

	if (criterion=="AIC" & min(n/k)<=40) {
		message( " [!] n/K<=40 for at least one model, consider using AICc criterion")
	}

	best.sum <- summary(pst.best)
	
	message(" [>] model ", id.best, " selected : ", criterion, "=", round(IC[id.best],2) , 
		" (C=", round(C[id.best],2), ")")
	message(" [>] ", best.sum@nodes, " nodes, ", best.sum@leaves, " leaves, ", 
		best.sum@freepar, " free parameters")
	
	if (output=="PST") {
		return(pst.best)
	} else if (output=="stats") {
		support <- rep(" ", nbmod)
		support[id.best] <- "***"
		tmp <- IC-IC[id.best]
		support[tmp>0 & tmp<=2] <- "**"
		support[tmp>2 & tmp<10] <- "*"

		res <- data.frame(Model=1:length(C), C=C, Nodes=nodes, Leaves=leaves, 
			Freepar=k, IC, Support=support)
		names(res)[6] <- criterion
		
		return(res)
	}
}
)
