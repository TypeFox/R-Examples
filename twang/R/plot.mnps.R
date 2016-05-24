##require(twang); data(AOD); AOD$crimjust[198:202] <- NA; mnps.AOD <- mnps(treat ~ illact + crimjust + subprob + subdep + white, data = AOD, estimand = "ATE", stop.method = c("ks.max","es.max"), n.trees = 1000, treatATT = 'community')
plot.mnps <- function(x,plots="optimize", pairwiseMax = TRUE, figureRows = NULL, color = TRUE, subset = NULL, treatments = NULL, singlePlot = NULL, multiPage = FALSE, time = NULL, print = TRUE, ...){
	
	stop.method <- tmt1 <- tmt2 <- sig <- NULL   

	if(!is.numeric(subset) & !is.null(subset)){
		if(!all(subset %in% x$stopMethods)) stop("The \"subset\" arugment must be NULL, numeric, or one of the stop.methods specified when fitting the mnps object.")
	}
   	if(is.null(subset)) subset <- 1:length(x$stopMethods)
   	if(!is.numeric(subset)){
   		hldLen <- 1:length(x$stopMethods)
   		subset <- hldLen[x$stopMethods %in% subset]
   	}
   	
   	subset <- sort(subset)
   	
   	if(is.null(subset)) stop("The \"subset\" arugment must either be NULL, numeric, or some subset of", print(x$stopMethods)) 
   	
   	if(length(treatments) > 2 & x$estimand == "ATE") stop("The \"treatments\" argument must be null or have length 1 or 2.")   	
   	
   	if(plots == 3 | plots == 4 | plots == 5 | plots == "es" | plots == "t" | plots == "ks"){
	   	if(pairwiseMax & !is.null(singlePlot)) warning("For this figure \"singlePlot\" argument is ignored when pairwiseMax = TRUE.")
	   	}
   	
   	if((length(treatments) > 1) & x$estimand == "ATT"){
   		warning("treatments argument must be null or have length 1 when estimand = ATT")
   	}
   	
   	if(pairwiseMax & !is.null(treatments)){
   		warning("treatments argument is ignored when pairwiseMax = TRUE.")
   	}
      
   	ltBl <- ifelse(color, "lightblue","gray80")
	rdCol <- ifelse(color, "red","black")
	stripBgCol <- ifelse(color, "#ffe5cc", "transparent")
	ptSymCol <- ifelse(color, "#0080ff", "black")
   
   ifelse(x$estimand == "ATE", noKS <- TRUE, noKS <- FALSE)
   
   subst <- whichVar <- pVal <- weighted <- NULL 
   
   if(length(plots) > 1) stop("The `plots' argument must be of length 1.")
   if(!(plots %in% c(1,2,3,4,5)) & !(plots %in% c("boxplot","ks","optimize","es","t")))
   stop("Invalid choice of `plots' argument.")
   
   if(plots == 2 | plots == "boxplot"){
   	
   	boxplot(x, color = color, stop.method = subset, multiPage = multiPage, singlePlot = singlePlot, time = time, print = print, ...)
   	
   }
   
   else if(plots == "optimize" | plots == 1){
   nPlot <- x$nFits
   
   ptHld <- vector(mode = "list", length = nPlot)
   for(i in 1:nPlot){
   	if(x$estimand == "ATT") ptNm <- paste("Balance for", x$levExceptTreatATT[i], "versus unweighted", x$treatATT)
   	else ptNm <- paste("Balance for", x$treatLev[i], "against others")
   	if(!is.null(time)){ptNm <- paste(ptNm, " (time ", time, ")", sep = "")}
   	ptHld[[i]] <- plot(x$psList[[i]], main = ptNm, plots = plots, noKS = TRUE, color = color, subset=subset)
   	}
   
   if(print) displayPlots(ptHld, figureRows = figureRows, singlePlot = singlePlot, multiPage = multiPage)
   return(ptHld)
   	
   }
   
   else if(!pairwiseMax){
   	if(x$estimand == "ATE"){
   			pairs <- bal.table(x)
	stp <- x$stopMethods
	unwTab <- subset(pairs, stop.method == "unw")
	hdUnwTab <- unwTab
	if(length(stp) > 1){
		for(i in 1:(length(stp) - 1)) unwTab <- rbind(unwTab, hdUnwTab)
	}
	
   	if(length(treatments) == 2) nPlotsTot <- 1
   	else if(length(treatments) == 1) nPlotsTot <- length(x$treatLev) - 1
   	else if(is.null(treatments)) nPlotsTot <- choose(length(x$treatLev), 2)
	
	wghtTab <- subset(pairs, stop.method != "unw")
	
	unwTab <- data.frame(unwTab, whichVar = 1:nrow(unwTab))
	wghtTab <- data.frame(wghtTab, whichVar = 1:nrow(wghtTab))
	
	plotTab <- rbind(unwTab, wghtTab)
	
	plotTab$stopMeth <- wghtTab$stop.method
	
	plotTab$weighted <- factor(rep(c("Unweighted","Weighted"), each = (nrow(plotTab)/2)))
	
	plotTab <- subset(plotTab, as.factor(plotTab$stopMeth) %in% levels(as.factor(plotTab$stopMeth))[subset])
	
	
	if(plots == 3 | plots == "es"){
	
	plotTab$bigger <- rep(as.numeric(wghtTab$std.eff.sz > unwTab$std.eff.sz),2)
	
	plotTab$sig <- as.numeric(plotTab$p < 0.05)
	
#   	if(is.null(subset))
#   	subset <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	yMax <- min(3,max(plotTab$std.eff.sz, na.rm=TRUE)) + .05	
	
   	nullPlot <- TRUE
   	
   	if(max(plotTab$std.eff.sz, na.rm=TRUE) > 3)
   	warning("Some effect sizes are larger than 3 and may not have been plotted.\n")	
   	
   	}
   	
   	allDat <- plotTab
   	if(!is.null(treatments) & !(all(treatments %in% x$treatLev))) {
   		print(x$treatLev)
   		stop("All elements of the \"treatments\" argument must be levels of the treatment variable, as printed above.")
   		}
   	if(length(treatments) == 1) plotTab <- subset(plotTab, tmt1 == treatments | tmt2 == treatments)
   	if(length(treatments) == 2) plotTab <- subset(plotTab, (tmt1 %in% treatments) & (tmt2 %in% treatments))	
   	

   	
   	ptNames <- ptList <- vector(mode = "list", length = nPlotsTot)
   	if(is.null(treatments)){
   		cnt <- 1
   		for(i in 1:(length(x$treatLev)-1))
   			for(j in (i+1):length(x$treatLev)){
   				ptNames[[cnt]] <- c(x$treatLev[i], x$treatLev[j]); cnt <- cnt + 1
   				}
   	}
   	else if(length(treatments) == 1){
   		subTreat <- x$treatLev[x$treatLev != treatments]
   		for(i in 1:nPlotsTot) ptNames[[i]] <- c(subTreat[i], treatments)
   	}
   	else if(length(treatments) == 2) ptNames[[1]] <- treatments
   	
   	
   	if(plots == 3 | plots == "es"){
   	
   	superPlotTab <- plotTab
   	
   	for(i in 1:nPlotsTot){
   		if(nPlotsTot > 1) plotTab <- subset(superPlotTab, tmt1 %in% ptNames[[i]] & tmt2 %in% ptNames[[i]])
	   	subsetHold <- ! plotTab$bigger 
   	
   		if(any(subsetHold)){   	
   			esDatTmp <- plotTab
   			esDatTmp$std.eff.sz[!subsetHold] <- NA	
   			ptNm <- paste("Balance of ", ptNames[[i]][1], " versus ", ptNames[[i]][2], sep = "")
   			if(!is.null(time)) ptNm <- paste(ptNm, " (time ", time, ")", sep = "")
   			pt1.1 <- xyplot(std.eff.sz ~ weighted | stopMeth, groups = whichVar, data = esDatTmp, 
   				scales = list(alternating = 1), ylim = c(-.05, yMax), type = "l", col = ltBl, 
   				as.table = TRUE, subset = subsetHold, par.settings = list(strip.background = list(col=stripBgCol)),
   				ylab = "Absolute standard difference", xlab = NULL, 
   				main = ptNm,
   				panel = function(...){
   					panel.abline(h=c(.2,.5,.8), col="gray80")
   					panel.xyplot(...)
		   		})
   			ptHold <- pt1.1
   			nullPlot <- FALSE
   		}
   	
	   	subsetHold <- plotTab$bigger & (as.factor(plotTab$stopMeth) %in% levels(as.factor(plotTab$stopMeth))[subset])
#	   	subsetHold <- plotTab$bigger & (as.factor(plotTab$stopMeth) %in% levels(as.factor(plotTab$stopMeth))[subset])

   	
	   	if(any(subsetHold)){
   		   	esDatTmp <- plotTab
   			esDatTmp$std.eff.sz[!subsetHold] <- NA		
   			ptNm <- paste("Balance of ", ptNames[[i]][1], " versus ", ptNames[[i]][2], sep = "")
   			if(!is.null(time)) ptNm <- paste(ptNm, " (time ", time, ")", sep = "")
   			pt1.2 <- xyplot(std.eff.sz ~ weighted | stopMeth, groups = whichVar, 
   			main = ptNm,
   			data = esDatTmp, ylab = "Absolute standard difference", xlab = NULL, as.table = TRUE, 
   			ylim = c(-.05, yMax), type = "l", col = rdCol, par.settings = list(strip.background = list(col=stripBgCol)),
   			lwd = 2)
   			if(nullPlot){
   				ptHold <- pt1.2
   				nullPlot <- FALSE
   			}
   			else {
   				ptHold <- ptHold + as.layer(pt1.2)
   				}
   		}
   	
#	   	subsetHold <- as.factor(plotTab$stopMeth) %in% levels(as.factor(plotTab$stopMeth))[subset]

   	
   		if(all(plotTab$p < 0.05, na.rm=TRUE)) pchHold <- 19
   		else if(all(plotTab$p >= 0.05, na.rm=TRUE)) pchHold <- 1
   		else pchHold <- c(1,19)
   	
		ptNm <- paste("Balance of ", ptNames[[i]][1], " versus ", ptNames[[i]][2], sep = "")
   		if(!is.null(time)) ptNm <- paste(ptNm, " (time ", time, ")", sep = "")

   		pt2 <- xyplot(std.eff.sz ~ weighted | stopMeth, groups = sig >= .05, data = plotTab,
   			ylab = "Absolute standard difference", xlab = NULL, 
   			ylim = c(-.05, yMax), type = "p", col = rdCol, pch = pchHold,
   			main = ptNm,
   			#subset = subsetHold, 
   			par.settings = list(strip.background = list(col=stripBgCol)))
   		ptHold <- ptHold + pt2
   	
   		pt1 <- ptHold
   		
   		ptList[[i]] <- pt1
   		
   		}
   		
   		}
   		
   	if(plots == 4 | plots == "t"){
  	nVar <- length(unique(allDat$var))
   	for(i in 1:nPlotsTot){
   		plotTab <- subset(allDat, tmt1 %in% ptNames[[i]] & tmt2 %in% ptNames[[i]])
   		plotTab$tRank <- NA
   	cnt <- 0
   	for(k in subset){
   		cnt <- cnt + 1
   	for(j in 1:length(x$psList)){	
		collapsed <- plotTab$p[plotTab$stopMeth == x$stopMethods[k] & plotTab$weighted == "Weighted"]	
		#collapsedUnw <- plotTab$p[plotTab$stopMeth == x$stopMethods[k] & plotTab$weighted == "Unweighted"]		
		collRanks <- rank(collapsed, ties.method = "first")
		#esBigHold <- collapsed > collapsedUnw
		#plotTab$tPVal[(1:(2*nVar)) + (cnt - 1)*2*nVar] <- c(collapsed, collapsedUnw)
		plotTab$tRank[plotTab$stopMeth == x$stopMethods[k] & plotTab$weighted == "Weighted"] <- collRanks
		plotTab$tRank[plotTab$stopMeth == x$stopMethods[k] & plotTab$weighted == "Unweighted"] <- collRanks

	}
	}
   	
   	   #	if(is.null(subst))	subst <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	n.var2 <- max(plotTab$tRank * (!is.na(plotTab$p)), na.rm=TRUE)
   	
	ptNm <- paste("Comparison of ", ptNames[[i]][1], " and ", ptNames[[i]][2], sep = "")
   	if(!is.null(time)) ptNm <- paste(ptNm, " (time ", time, ")", sep = "")
   	
   	pt1 <- xyplot(p~tRank|stopMeth, groups = weighted, data=plotTab, xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "t- and chi-squared p-values", pch = c(19,1), col = "black", scales = list(alternating = 1), par.settings = list(strip.background = list(col=stripBgCol)),
   	main = ptNm,
   	ylim = c(-.1, 1.1), 
   	   	panel = function(...){
   	   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col=ltBl, type="l")
#   		panel.abline(a=0, b=1, col="lightblue")
   		panel.xyplot(...)
   		}
)

	ptList[[i]] <- pt1
   		
   	}
   	
   	}
   	
   	if(plots == 5 | plots == "ks"){

  	nVar <- length(unique(allDat$var))
   	for(i in 1:nPlotsTot){
   		plotTab <- subset(allDat, tmt1 %in% ptNames[[i]] & tmt2 %in% ptNames[[i]])
   		plotTab$ksRank <- NA
   	cnt <- 0
   	for(k in subset){
   		cnt <- cnt + 1
   	for(j in 1:length(x$psList)){	
		collapsed <- plotTab$ks.pval[plotTab$stopMeth == x$stopMethods[k] & plotTab$weighted == "Weighted"]	
		collRanks <- rank(collapsed, ties.method = "first")
		plotTab$ksRank[plotTab$stopMeth == x$stopMethods[k] & plotTab$weighted == "Weighted"] <- collRanks
		plotTab$ksRank[plotTab$stopMeth == x$stopMethods[k] & plotTab$weighted == "Unweighted"] <- collRanks

	}
	}
   	
   	   #	if(is.null(subst))	subst <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	n.var2 <- max(plotTab$ksRank * (!is.na(plotTab$ks.pval)), na.rm=TRUE)
   	
	ptNm <- paste("Comparison of ", ptNames[[i]][1], " and ", ptNames[[i]][2], sep = "")
   	if(!is.null(time)) ptNm <- paste(ptNm, " (time ", time, ")", sep = "")
   	
   	
   	pt1 <- xyplot(ks.pval~ksRank|stopMeth, groups = weighted, data=plotTab, xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "KS test p-values", pch = c(19,1), col = "black", scales = list(alternating = 1), par.settings = list(strip.background = list(col=stripBgCol)),
   	main = ptNm,
   	ylim = c(-.1, 1.1),
   	   	panel = function(...){
   	   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col=ltBl, type="l")
#   		panel.abline(a=0, b=1, col="lightblue")
   		panel.xyplot(...)
   		}
)

	ptList[[i]] <- pt1
   		
   	}
   		
   		
   	}	
	
	if(nPlotsTot == 1) return(pt1)
	else {
		if(print) displayPlots(ptList, figureRows = figureRows, singlePlot = singlePlot, multiPage = multiPage)	
		return(ptList)
	}
	

   		
   	}
   else{
   nPlot <- x$nFits
   ptHld <- vector(mode = "list", length = nPlot)
   for(i in 1:nPlot){
   	if(x$estimand == "ATT") ptNm <- paste("Balance for", x$levExceptTreatATT[i], "versus unweighted", x$treatATT)
   	else ptNm <- paste("Balance for", x$treatLev[i], "against others")
   	if(!is.null(time)) ptNm <- paste(ptNm, " (time ", time, ")", sep = "")   	
   	ptHld[[i]] <- plot(x$psList[[i]], main = ptNm, plots = plots, color = color, subset=subset)
   }
   
	if(!is.null(treatments)){
   		ptNum <- 1:nPlot[x$levExceptTreatAtt == treatments]
   		return(ptHld[[ptNum]])
   	}
   	else {
   		if(print) displayPlots(ptHld, figureRows = figureRows, singlePlot = singlePlot, multiPage = multiPage)
   		return(ptHld)
   	}   
   

}
}

   else{  ## if pairwiseMax and plots == something other than optimize
   	
   	n.tp <- length(x$psList[[1]]$desc)
   	n.psFits <- length(x$psList)
   	
   if (plots == "es" || plots == 3)	{ ## es plot
   	
#   	esDat <- makePlotDat(x$psList[[1]], whichPlot = 3, subsetStopMeth = subset)
#   	nVar <- length(makePlotDat(x$psList[[1]], whichPlot = 3, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE)$pVal)
#   	effSzList <- effSzListUnw <- pValListUnw <- pValList <- matrix(NA, nrow = nVar, ncol = length(x$psList))\
   	   	
   	pwc <- bal.table(x, collapse.to = "covariate")
   	
   	nVar <- nrow(subset(pwc, stop.method == "unw"))
   	
   	hldEffSz <- hldPVal <- hldEsBig <- hldWhichComp <- hldWhichVar <- hldWeighted <- rep(NA, (nVar * length(subset) * 2))

   	
   	cnt <- 0
   	for(k in subset){
   		cnt <- cnt + 1
			
		collapsed <- pwc$max.std.eff.sz[pwc$stop.method == x$stopMethods[k]]
   		collapsedP <- pwc$min.p[pwc$stop.method == x$stopMethods[k]] < .05	
   		collapsedUnw <- pwc$max.std.eff.sz[pwc$stop.method == "unw"]
   		collapsedUnwP <- pwc$min.p[pwc$stop.method == "unw"] < .05		
	
		esBigHold <- collapsed > collapsedUnw
#		esDat$effectSize[(1:(2*nVar)) + (cnt-1)*2*nVar] <- c(collapsed, collapsedUnw)
#		esDat$pVal[(1:(2*nVar)) + (cnt-1)*2*nVar] <- c(collapsedP, collapsedUnwP)
#		esDat$esBig[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(esBigHold, 2)
		hldEffSz[(1:(2*nVar)) + (cnt-1)*2*nVar] <- c(collapsed, collapsedUnw)
		hldPVal[(1:(2*nVar)) + (cnt-1)*2*nVar] <- c(collapsedP, collapsedUnwP)
		hldEsBig[(1:(2*nVar)) + (cnt-1)*2*nVar] <-  rep(esBigHold, 2)
		hldWhichComp[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(x$stopMethods[k], 2 * nVar)
		hldWhichVar[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(1:nVar, 2)
		hldWeighted[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(c("Weighted","Unweighted"), each = nVar)
	}
   	
   	esDat <- data.frame(effectSize = abs(hldEffSz), esBig = hldEsBig, whichComp = hldWhichComp, weighted = hldWeighted, whichVar = hldWhichVar, pVal = hldPVal)
   	
   	yMax <- min(3,max(esDat[,1])) + .05	
   	
   	if(max(esDat[,1], na.rm=TRUE) > 3)
   	warning("Some effect sizes are larger than 3 and may not have been plotted.\n")	
   	
    
    nullPlot <- TRUE
       	
   	subsetHold <- !esDat$esBig 
   	
   	ptNm <- NULL
   	if(!is.null(time)) ptNm <- paste("Time ", time, sep = "")
   	
   	if(any(subsetHold)){
   		esDatTmp <- esDat
   		esDatTmp$effectSize[!subsetHold] <- NA
   		pt1.1 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, data = esDatTmp, scales = list(alternating = 1),
   			ylim = c(-.05, yMax), type = "l", col = ltBl, as.table = TRUE, main = ptNm, 
   			ylab = "Absolute standardized difference \n (maximum pairwise)", xlab = NULL, par.settings = list(strip.background = list(col=stripBgCol)),
   			panel = function(...){
   				panel.abline(h=c(.2,.5,.8), col="gray80")
   				panel.xyplot(...)
		   	})
   		nullPlot <- FALSE
   		currPt <- pt1.1
   	}
   	
   	subsetHold <- esDat$esBig 
   	
   	if(any(subsetHold)){
   		esDatTmp <- esDat
   		esDatTmp$effectSize[!subsetHold] <- NA
   		pt1.2 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, 
   			data = esDatTmp, ylab = "Absolute standard difference", xlab = NULL, as.table = TRUE,
   			ylim = c(-.05, yMax), type = "l", col = rdCol, par.settings = list(strip.background = list(col=stripBgCol)), main = ptNm,
   			lwd = 2)
   		if(!nullPlot){
   			currPt <- currPt + pt1.2
   		}
   		else{
   			currPt <- pt1.2
   			nullPlot <- FALSE
   		}
   	}
   	
   	
   	if(all(esDat$pVal >= .05)) pchHold <- 19
   	else if(all(esDat$pVal < .05)) pchHold <- 1
   	else pchHold <- c(19,1)
   	
   	pt2 <- xyplot(effectSize ~ weighted | whichComp, groups = (pVal < 0.05), data = esDat,
   	ylab = "Absolute standard difference", xlab = NULL, as.table = TRUE, main = ptNm, 
   	ylim = c(-.05, yMax), type = "p", col = rdCol, pch = pchHold,par.settings = list(strip.background = list(col=stripBgCol)))
   	if(!nullPlot) currPt <- currPt + pt2
   	else currPt <- pt2
   	
   	return(currPt)
   	
   	} 			
   
   if (plots == "t" || plots == 4) { ## t p-values plot

#   	esDat <- makePlotDat(x$psList[[1]], whichPlot = 4, subsetStopMeth = subset)
#   	nVar <- length(makePlotDat(x$psList[[1]], whichPlot = 4, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE))
#   	effSzList <- effSzListUnw <- matrix(NA, nrow = nVar, ncol = length(x$psList))   	
   	pwc <- bal.table(x, collapse.to = "covariate")
   	
   	nVar <- nrow(subset(pwc, stop.method == "unw"))
   	
   	hldTPVal <- hldTRank <- hldWhichComp <- hldWhichVar <- hldWeighted <- rep(NA, (nVar * length(subset) * 2))
   	
   	ptNm <- NULL
   	if(!is.null(time)) ptNm <- paste("Time ", time, sep = "")
  	
   	
   	cnt <- 0
   	for(k in subset){
   		cnt <- cnt + 1
   	for(j in 1:length(x$psList)){	
		collapsed <- pwc$min.p[pwc$stop.method == x$stopMethods[k]]	
		collapsedUnw <- pwc$min.p[pwc$stop.method == 'unw']		
		collRanks <- rank(collapsed, ties.method = "first")
#		esBigHold <- collapsed > collapsedUnw
		hldTPVal[(1:(2*nVar)) + (cnt - 1)*2*nVar] <- c(collapsed, collapsedUnw)
		hldTRank[(1:(2*nVar)) + (cnt - 1)*2*nVar] <- rep(collRanks, 2)
		hldWhichComp[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(x$stopMethods[k], 2 * nVar)
		hldWhichVar[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(1:nVar, 2)
		hldWeighted[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(c("Weighted","Unweighted"), each = nVar)		
	}
	}
	
	esDat <- data.frame(tPVal = hldTPVal, tRank = hldTRank, whichComp = hldWhichComp, weighted = hldWeighted)
   	   	
   	n.var2 <- max(esDat$tRank * (!is.na(esDat$tPVal)), na.rm=TRUE)
   	
   	pt1 <- xyplot(tPVal~tRank|whichComp, groups = weighted, data=esDat, xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "t- and chi-squared p-values \n (pairwise minimum)", pch = c(19,1), col = "black", scales = list(alternating = 1), par.settings = list(strip.background = list(col=stripBgCol)), main = ptNm,
   	#subst = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subst]) & (esDat$tRank <= n.var2), 
   	ylim = c(-.1, 1.1), 
   	   	panel = function(...){
   	   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col=ltBl, type="l")
#   		panel.abline(a=0, b=1, col="lightblue")
   		panel.xyplot(...)
   		}
)
   	
   	}
   
   if (plots =="ks" || plots ==5) {  ## ks plot
#   	if(x$estimand =="ATE") pwc <- pairwiseComparison(x, collapse.to = "covariate")

   	pwc <- bal.table(x, collapse.to = "covariate")
   	
   	nVar <- nrow(subset(pwc, stop.method == "unw"))

   	ptNm <- NULL
   	if(!is.null(time)) ptNm <- paste("Time ", time, sep = "")
   	


#   	esDat <- makePlotDat(x$psList[[1]], whichPlot = 5, subsetStopMeth = subset)
#   	nVar <- length(makePlotDat(x$psList[[1]], whichPlot = 5, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE))
#   	effSzList <- effSzListUnw <- matrix(NA, nrow = nVar, ncol = length(x$psList))  

#   	nVar <- nrow(subset(pwc, stop.method == "unw"))
   	
   	hldKsPVal <- hldKsRank <- hldWhichComp <- hldWhichVar <- hldWeighted <- rep(NA, (nVar * length(subset) * 2))
    	
   	
   	cnt <- 0
   	for(k in subset){
   		cnt <- cnt + 1
   	for(j in 1:length(x$psList)){	
#	x2 <- x$psList[[j]]
#	effSzList[,j] <- makePlotDat(x2, whichPlot = 5, subsetStopMeth = k, yOnly = TRUE, incUnw = FALSE)
#	effSzListUnw[,j] <- makePlotDat(x2, whichPlot = 5, subsetStopMeth = 1, yOnly = TRUE, incUnw = TRUE)[1:nVar]
#	collapsed <- apply(effSzList, 1, max)
#	if(x$estimand == "ATE") collapsed <- pwc$min.ks.pval[pwc$stop.method == x$stopMethods[k]]
#	collapsedUnw <- apply(effSzListUnw, 1, max)
	collapsed <- pwc$min.ks.pval[pwc$stop.method == x$stopMethods[k]]	
	collapsedUnw <- pwc$min.ks.pval[pwc$stop.method == 'unw']			
#	if(x$estimand == "ATE") collapsedUnw <- pwc$min.ks.pval[pwc$stop.method == "unw"]	
	collRanks <- rank(collapsed, ties.method = "first")
#	esBigHold <- collapsed > collapsedUnw
		hldKsPVal[(1:(2*nVar)) + (cnt - 1)*2*nVar] <- c(collapsed, collapsedUnw)
		hldKsRank[(1:(2*nVar)) + (cnt - 1)*2*nVar] <- rep(collRanks, 2)
		hldWhichComp[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(x$stopMethods[k], 2 * nVar)
		hldWhichVar[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(1:nVar, 2)
		hldWeighted[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(c("Weighted","Unweighted"), each = nVar)		}
	} 

	esDat <- data.frame(ksPVal = hldKsPVal, ksRank = hldKsRank, whichComp = hldWhichComp, weighted = hldWeighted)

   	
   	if(is.null(subst))	subst <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	n.var2 <- max(esDat$ksRank*(!is.na(esDat$ksPVal)), na.rm=TRUE)
   	
   	pt1 <- xyplot(ksPVal~ksRank|whichComp, groups=weighted, scales = list(alternating = 1), data = esDat,ylim = c(-.1, 1.1),
   	 xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "KS test p-values", pch = c(19,1), col="black", par.settings = list(strip.background = list(col=stripBgCol)), main = ptNm, 
   	subset = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subst]) & (esDat$ksRank <= n.var2),
   	panel = function(...){
   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col=ltBl, type="l")
   		panel.xyplot(...)
   	})
   	}
   	
return(pt1)
	
	}

}

