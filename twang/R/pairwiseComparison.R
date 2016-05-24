pairwiseComparison <- function(x, collapse.to = c("pair","covariate","stop.method")[1], na.action = "level"){
	if(class(x) != "mnps") stop("pairwiseComparison only defined for mnps objects fit with estimand = \"ATE\"")
	if(x$estimand != "ATE") stop("pairwiseComparison only defined for mnps objects fit with estimand = \"ATE\"")	
	stop.method <- NULL
	treatLevs <- x$treatLev
	nTreat <- length(treatLevs)
	treatInds1 <- NULL
	for(i in 1:(nTreat-1)) treatInds1 <- c(treatInds1, rep(i,nTreat - i))
	treatInds2 <- NULL
	for(i in 2:(nTreat)) treatInds2 <- c(treatInds2, i:nTreat)
	subDt <- x$data
	subDt[,x$treat.var] <- as.numeric(subDt[,x$treat.var] == levels(subDt[,x$treat.var])[1])	
	tabForNames <- desc.wts(subDt, w = rep(1,nrow(subDt)), sampw = rep(1,nrow(subDt)), vars = x$psList[[1]]$gbm.obj$var.names, treat.var = x$treat.var, na.action = "level", verbose = FALSE, alerts.stack = 0, estimand = x$estimand, multinom = FALSE, fillNAs = TRUE)$bal.tab$results
	nRowBalTab <- nrow(tabForNames)
	rwNms <- row.names(tabForNames)
	
	hldBalTabs <- vector(mode = "list", length = length(treatInds1))
	for(i in 1:length(treatInds1)){
		subDt <- x$data[x$data[,x$treat.var] %in% treatLevs[c(treatInds1[i], treatInds2[i])], ]
		subDt[,x$treat.var] <- as.numeric(subDt[,x$treat.var] == treatLevs[treatInds1[i]])
		bTab <- desc.wts(subDt, w = rep(1,nrow(subDt)), sampw = rep(1,nrow(subDt)), vars = x$psList[[1]]$gbm.obj$var.names, treat.var = x$treat.var, na.action = "level", verbose = FALSE, alerts.stack = 0, estimand = x$estimand, multinom = FALSE, fillNAs = TRUE)$bal.tab$results
		bTab <- bTab[rwNms,]
		bTab[is.na(bTab$tx.mn),] <- ifelse(names(bTab) %in% c("p","ks.pval"), 1, 0)
		hldBalTabs[[i]] <- bTab
	}
	
	allTabs <- do.call(rbind, hldBalTabs)
	
	
	
	stpMth <- c(x$stopMethods)
	
	for(i in 1:length(stpMth)){
		hldBalTabs <- vector(mode = "list", length = length(treatInds1))
		wgt <- get.weights(x, stop.method = stpMth[i], estimand = x$estimand)
		for(i in 1:length(treatInds1)){
			subDt <- x$data[x$data[,x$treat.var] %in% treatLevs[c(treatInds1[i], treatInds2[i])], ]
			subW <- wgt[x$data[,x$treat.var] %in% treatLevs[c(treatInds1[i], treatInds2[i])]]
			subDt[,x$treat.var] <- as.numeric(subDt[,x$treat.var] == treatLevs[treatInds1[i]])
			bTab <- desc.wts(subDt, w = subW, sampw = rep(1,nrow(subDt)), vars = x$psList[[1]]$gbm.obj$var.names, treat.var = x$treat.var, na.action = "level", verbose = FALSE, alerts.stack = 0, estimand = x$estimand, multinom = FALSE, fillNAs = TRUE)$bal.tab$results
			bTab <- bTab[rwNms,]
			bTab[is.na(bTab$tx.mn),] <- ifelse(names(bTab) %in% c("p","ks.pval"), 1, 0)
			hldBalTabs[[i]] <- bTab
		}
	
		allTabs <- rbind(allTabs, rbind(do.call(rbind, hldBalTabs)))
		
	}
	
	btb <- bal.table.hidden(x, pairwise = FALSE, digits = 10)
	balTab1 <- btb[[1]][[1]]
	
	tableShell <- data.frame(tmt1 = rep(rep(treatLevs[treatInds1], each = nRowBalTab), length(stpMth) + 1), tmt2 = rep(rep(treatLevs[treatInds2], each = nRowBalTab), length(stpMth) + 1), var = rep(rwNms, length(treatInds1) * (length(stpMth) + 1)), mean1 = NA, mean2 = NA, pop.sd = rep(balTab1$ct.sd, length(treatInds1) * (length(stpMth) + 1)), std.eff.sz = NA, p = NA, ks = NA, ks.pval = NA, stop.method = rep(c("unw", stpMth), each = (length(treatInds1) * length(rwNms))))
	
	
#	btb <- bal.table(x)
#	hldMn1 <- hldMn2 <- NULL
#	for(j in 1: length(treatInds1)){
#		hldMn1 <- c(hldMn1, btb[[treatInds1[j]]][["unw"]]$tx.mn)
#		hldMn2 <- c(hldMn2, btb[[treatInds2[j]]][["unw"]]$tx.mn)			
#	}	
#	for(i in 1:length(stpMth)){
#		for(j in 1: length(treatInds1)){
#			hldMn1 <- c(hldMn1, btb[[treatInds1[j]]][[stpMth[i]]]$tx.mn)
#			hldMn2 <- c(hldMn2, btb[[treatInds2[j]]][[stpMth[i]]]$tx.mn)			
#		}
#	}
	
	tableShell$mean1 <- allTabs$tx.mn; tableShell$mean2 <- allTabs$ct.mn
	
	tableShell$std.eff.sz <- abs(tableShell$mean1 - tableShell$mean2)/tableShell$pop.sd
	
	tableShell$p <- allTabs$p; tableShell$ks <- allTabs$ks; tableShell$ks.pval <- allTabs$ks.pval
	
	if(collapse.to == "pair") return(tableShell)
	else {
		subDt <- subset(tableShell, stop.method == "unw")
		asmds <- apply(matrix(subDt$std.eff.sz, ncol = length(treatInds1)), 1, max, na.rm = TRUE)
		pvals <- apply(matrix(subDt$p, ncol = length(treatInds1)), 1, min, na.rm = TRUE)
		kss <- apply(matrix(subDt$ks, ncol = length(treatInds1)), 1, max, na.rm = TRUE)		
		kspvals <- apply(matrix(subDt$ks.pval, ncol = length(treatInds1)), 1, min, na.rm = TRUE)		
		for(i in 1:length(stpMth)){
			subDt <- subset(tableShell, stop.method == stpMth[i])
			asmds <- c(asmds, apply(matrix(subDt$std.eff.sz, ncol = length(treatInds1)), 1, max, na.rm = TRUE))
			pvals <- c(pvals, apply(matrix(subDt$p, ncol = length(treatInds1)), 1, min, na.rm = TRUE))
			kss <- c(kss, apply(matrix(subDt$ks, ncol = length(treatInds1)), 1, max, na.rm = TRUE))		
			kspvals <- c(kspvals, apply(matrix(subDt$ks.pval, ncol = length(treatInds1)), 1, min, na.rm = TRUE))		

		}

		redTableShell <- data.frame(var = rep(rwNms, length(stpMth) + 1), max.std.eff.sz = asmds, min.p = pvals, max.ks = kss, min.ks.pval = kspvals, stop.method = rep(c("unw", stpMth), each = length(rwNms)))
		if(collapse.to == "covariate") return(redTableShell)
		else {
			subDt <- subset(redTableShell, stop.method == "unw")
			asmds <- max(subDt$max.std.eff.sz, na.rm = TRUE)
			pvals <- min(subDt$min.p, na.rm = TRUE)
			kss <- max(subDt$max.ks, na.rm = TRUE)
			kspvals <- min(subDt$min.ks.pval, na.rm = TRUE)
	
		
			for(i in 1:length(stpMth)){
				subDt <- subset(redTableShell, stop.method == stpMth[i])
				asmds <- c(asmds, max(subDt$max.std.eff.sz, na.rm = TRUE))
				pvals <- c(pvals, min(subDt$min.p, na.rm = TRUE))
				kss <- c(kss, max(subDt$max.ks, na.rm = TRUE))
				kspvals <- c(kspvals, min(subDt$min.ks.pval, na.rm = TRUE))
			}
		
			redRedTableShell <- data.frame(max.std.eff.sz = asmds, min.p = pvals, max.ks = kss, min.ks.pval = kspvals, stop.method = c("unw",stpMth))
			return(redRedTableShell)
		}
		
	}
	
}