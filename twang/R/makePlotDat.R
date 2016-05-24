makePlotDat <- function(x, whichPlot, subsetStopMeth=NULL, yOnly = FALSE, incUnw = TRUE){
	## whichPlot:  1 = optimize; 2 = boxplot -- Ignore for this fcn
	## 3 = esPlot; 4 = t p-values; 5 = ks p-values
	
#	n.tp <- ifelse(class(x) == "dxwts", length(x$desc), ncol(weights))
	if(class(x) == "mnps") {
		n.tp <- length(x$psList[[1]]$desc)
		n.var <- nrow(x$psList[[1]]$desc$unw$bal.tab$results)
		}

	else if(class(x) == "ps") {
		n.tp <- length(x$desc) 
		n.var <- nrow(x$desc$unw$bal.tab$results)
	}
	
	
	

	
	if(whichPlot == 1){
		
	longBal <- matrix(t(x$balance))
	optDat <- data.frame(balance = longBal, iteration = rep(x$iters, each = n.tp-1), stopRule = names(x$desc)[-1])

	return(optDat)
		 
		
		
	}
	else if(whichPlot == 3){
		tpHld <- NULL	
		
		desc.unw <- x$desc[[1]]
		hldESBig <- hldEffSz <- hldPVal <- pVal <- effectSize <- esBig <- NULL
		stopIndices <- 2:n.tp
		if(!is.null(subsetStopMeth)) stopIndices <- stopIndices[subsetStopMeth]
		for(i in stopIndices){
			desc.temp <- x$desc[[i]]
			iter <- desc.temp$n.trees
			tpHld <- c(tpHld, names(x$desc[i]))
			
			hldEffSz <- c(hldEffSz, desc.temp$bal.tab$results$std.eff.sz)
			hldESBig <- c(hldESBig, abs(desc.temp$bal.tab$results$std.eff.sz) > abs(desc.unw$bal.tab$results$std.eff.sz))
			hldPVal <- c(hldPVal, desc.temp$bal.tab$results$p)
			if(incUnw){
				hldEffSz <- c(hldEffSz, desc.unw$bal.tab$results$std.eff.sz)
				hldESBig <- c(hldESBig, abs(desc.temp$bal.tab$results$std.eff.sz) > abs(desc.unw$bal.tab$results$std.eff.sz))
				hldPVal <- c(hldPVal, desc.unw$bal.tab$results$p)
			}
			
		}
		
	   	whichComp = as.factor(rep(tpHld, each = 2*n.var)) 
		weighted = as.factor(rep(rep(c("Weighted", "Unweighted"), each = n.var), length(stopIndices))) 
		whichVar = factor(rep(1:n.var, 2*length(stopIndices)))  		

		if(!yOnly){
			esDat <- data.frame(effectSize = abs(hldEffSz), esBig = hldESBig, whichComp = whichComp, 
			weighted = weighted, whichVar = whichVar, pVal = hldPVal)
		}
		else{
		 esDat <- list(effectSize = hldEffSz, pVal = hldPVal)
		 }
		 
		 if(any(is.na(esDat$pVal))) 
		 	for(i in 2:length(esDat$pVal)) 
		 		if(is.na(esDat$pVal[i])) esDat$pVal[i] <- esDat$pVal[i-1]
		
		return(esDat)
		
	}
	else if(whichPlot == 4){
		tpHld <- NULL	

		desc.unw <- x$desc[[1]]
		hldTVal <- hldTRank <- tPVal <- tRank <- NULL
		stopIndices <- 2:n.tp		
		if(!is.null(subsetStopMeth)) stopIndices <- stopIndices[subsetStopMeth]
		for(i in stopIndices){
			desc.temp <- x$desc[[i]]
			iter <- desc.temp$n.trees
			tpHld <- c(tpHld, names(x$desc[i]))

			hldTVal <- c(hldTVal, desc.temp$bal.tab$results$p)
			hldTRank <- c(hldTRank, rank(desc.temp$bal.tab$results$p, ties.method = "first"))	
			if(incUnw){
				hldTVal <- c(hldTVal, desc.unw$bal.tab$results$p)
				hldTRank <- c(hldTRank, rank(desc.unw$bal.tab$results$p, ties.method = "first"))
			}
			
		whichComp = as.factor(rep(tpHld, each = 2*n.var)) 
		weighted = as.factor(rep(rep(c("Weighted", "Unweighted"), each = n.var), length(stopIndices))) 
		whichVar = factor(rep(1:n.var, 2*length(stopIndices)))  
					
		}
		
		if(!yOnly){
			esDat <- data.frame(tPVal = hldTVal, tRank = hldTRank, whichComp = whichComp, 
			weighted = weighted, whichVar = whichVar)
		}
		else esDat <- unlist(data.frame(tPVal = hldTVal))	
	}
	else if(whichPlot == 5){
		tpHld <- NULL	
		
		desc.unw <- x$desc[[1]]
		hldksPVal <- hldksPValRanks <- ksPVal <- ksRank <- NULL
		stopIndices <- 2:n.tp		
		if(!is.null(subsetStopMeth)) stopIndices <- stopIndices[subsetStopMeth]
		for(i in stopIndices){
			desc.temp <- x$desc[[i]]
			iter <- desc.temp$n.trees
			tpHld <- c(tpHld, names(x$desc[i]))

			hldksPVal <- c(hldksPVal, desc.temp$bal.tab$results$ks.pval)
			hldksPValRanks <- c(hldksPValRanks, rank(desc.temp$bal.tab$results$ks.pval, ties.method = "first"))
			if(incUnw){
				hldksPVal <- c(hldksPVal, desc.unw$bal.tab$results$ks.pval)
				hldksPValRanks <- c(hldksPValRanks, rank(desc.unw$bal.tab$results$ks.pval, ties.method = "first"))
			}
						
		}
		
	   	whichComp = as.factor(rep(tpHld, each = 2*n.var)) 
		weighted = as.factor(rep(rep(c("Weighted", "Unweighted"), each = n.var), length(stopIndices))) 
		whichVar = factor(rep(1:n.var, 2*length(stopIndices)))  
		
		
		if(!yOnly){
			esDat <- data.frame(ksPVal = hldksPVal, ksRank = hldksPValRanks, whichComp = whichComp, 
			weighted = weighted, whichVar = whichVar)
		}
		else esDat <- unlist(data.frame(ksPVal = hldksPVal))			
		
	}
	
	
	}