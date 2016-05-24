##require(twang); data(AOD); AOD$crimjust[198:202] <- NA; mnps.AOD <- mnps(treat ~ illact + crimjust + subprob + subdep + white, data = AOD, estimand = "ATT", stop.method = c("ks.max","es.max"), n.trees = 1000, treatATT = 'community')
plot.mnps <- function(x,plots="optimize", pairwiseMax = TRUE, figureRows = 1, color = TRUE, subset = NULL, ...)
{
   # Creates diag.plot plots and sends to current device
   # x:     ps object 
   # label: Label added to the plot titles

   # extract the propensity scores and weights from the ps object
   

   	if(is.null(subset)) subset <- 1:length(x$stopMethods)
      
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
   	
   	boxplot(x, color = color, stop.method = subset, ...)
   	
   }
   
   else{
   
   if(!pairwiseMax | plots == "optimize" | plots == 1){
   
   nPlot <- x$nFits
   ptHld <- vector(mode = "list", length = nPlot)
   for(i in 1:nPlot){
   	if(x$estimand == "ATT") ptNm <- paste("Balance for", x$levExceptTreatATT[i], "versus unweighted", x$treatATT)
   	else ptNm <- paste("Balance for", x$treatLev[i], "against others")
   	ptHld[[i]] <- plot(x$psList[[i]], main = ptNm, plots = plots, noKS = TRUE, color = color, subset=subset, ...)
   }


figCol <- ceiling(nPlot/figureRows)

if(dev.cur() == 1) dev.new()

curCol <- curRow <- 1

for(i in 1:(nPlot-1)){
	print(ptHld[[i]], split = c(curCol,curRow,nx = figCol,ny = figureRows), more = TRUE)
	if(curCol < figCol){
		curCol <- curCol + 1
	}
	else {
		curCol <- 1
		curRow <- curRow + 1
	}
}

print(ptHld[[nPlot]], split = c(curCol,curRow,nx = figCol,ny = figureRows), more = FALSE)

}

else{  ## if pairwiseMax and plots == something other than optimize

n.tp <- length(x$psList[[1]]$desc)
n.psFits <- length(x$psList)
   	
#   	return(esDat)   
   	
   if (plots == "es" || plots == 3)	{ ## es plot
   	
   	esDat <- makePlotDat(x$psList[[1]], whichPlot = 3, subsetStopMeth = subset)
   	nVar <- length(makePlotDat(x$psList[[1]], whichPlot = 3, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE)$pVal)
   	effSzList <- effSzListUnw <- pValListUnw <- pValList <- matrix(NA, nrow = nVar, ncol = length(x$psList))
   	
   	pwc <- bal.table(x, collapse.to = "covariate")
   	
   	#if(x$estimand == "ATE"){
   	#	pwc <- pairwiseComparison(x, collapse.to = "covariate")	
   	#	}
   		

   		
   	
   	cnt <- 0
   	for(k in subset){
   		cnt <- cnt + 1
   		#if(x$estimand == "ATT"){
   		#for(j in 1:length(x$psList)){	
			#x2 <- x$psList[[j]]
			#effSzList[,j] <- makePlotDat(x2, whichPlot = 3, subsetStopMeth = k, yOnly = TRUE, incUnw = FALSE)$effectSize
			#effSzListUnw[,j] <- makePlotDat(x2, whichPlot = 3, subsetStopMeth = 1, yOnly = TRUE, incUnw = TRUE)$effectSize[(nVar + 1):(2*nVar)]
			#pValListUnw[,j] <- makePlotDat(x2, whichPlot = 3, subsetStopMeth = 1, yOnly = TRUE, incUnw = TRUE)$pVal[(nVar + 1):(2*nVar)]
			#pValList[,j] <- makePlotDat(x2, whichPlot = 3, subsetStopMeth = k, yOnly = TRUE, incUnw = FALSE)$pVal	
			#}
			#effSzList <- abs(effSzList)
			#effSzListUnw <- abs(effSzListUnw)
			#collapsed <- apply(effSzList, 1, max)
			#collapsedUnw <- apply(effSzListUnw, 1, max)
			#collapsedP <- apply(pValList < .05, 1, max)
			#collapsedUnwP <- apply(pValListUnw < .05, 1, max)	
			
			collapsed <- pwc$max.std.eff.sz[pwc$stop.method == x$stopMethods[k]]
   			collapsedP <- pwc$min.p[pwc$stop.method == x$stopMethods[k]] < .05
#   			collapsedUnw <- rep(pwc$max.std.eff.sz[pwc$stop.method == "unw"], length(subset))
#   			collapsedUnwP <- rep(pwc$min.p[pwc$stop.method == "unw"], length(subset)) < .05		
   			collapsedUnw <- pwc$max.std.eff.sz[pwc$stop.method == "unw"]
   			collapsedUnwP <- pwc$min.p[pwc$stop.method == "unw"] < .05		
		#}
	
	
		esBigHold <- collapsed > collapsedUnw
		esDat$effectSize[(1:(2*nVar)) + (cnt-1)*2*nVar] <- c(collapsed, collapsedUnw)
		esDat$pVal[(1:(2*nVar)) + (cnt-1)*2*nVar] <- c(collapsedP, collapsedUnwP)
		esDat$esBig[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(esBigHold, 2)
	}
   	
   	
   	yMax <- min(3,max(esDat[,1])) + .05	
   	
   	if(max(esDat[,1], na.rm=TRUE) > 3)
   	warning("Some effect sizes are larger than 3 and may not have been plotted.\n")	
   	
    
    nullPlot <- TRUE
       	
   	subsetHold <- !esDat$esBig 
   	
   	
   	if(any(subsetHold)){
   		esDatTmp <- esDat
   		esDatTmp$effectSize[!subsetHold] <- NA
   		pt1.1 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, data = esDatTmp, scales = list(alternating = 1),
   			ylim = c(-.05, yMax), type = "l", col = ltBl, as.table = TRUE,
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
   			ylim = c(-.05, yMax), type = "l", col = rdCol, par.settings = list(strip.background = list(col=stripBgCol)),
   			lwd = 2, ...)
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
   	ylab = "Absolute standard difference", xlab = NULL, as.table = TRUE, 
   	ylim = c(-.05, yMax), type = "p", col = rdCol, pch = pchHold,par.settings = list(strip.background = list(col=stripBgCol)),
   	...)
   	if(!nullPlot) currPt <- currPt + pt2
   	else currPt <- pt2
   	
   	return(currPt)
   	
   	}
   						
 			
   
   if (plots == "t" || plots == 4) { ## t p-values plot

   	esDat <- makePlotDat(x$psList[[1]], whichPlot = 4, subsetStopMeth = subset)
   	nVar <- length(makePlotDat(x$psList[[1]], whichPlot = 4, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE))
   	effSzList <- effSzListUnw <- matrix(NA, nrow = nVar, ncol = length(x$psList))   	
   	pwc <- bal.table(x, collapse.to = "covariate")
  	
   	
   	cnt <- 0
   	for(k in subset){
   		cnt <- cnt + 1
   	for(j in 1:length(x$psList)){	
		collapsed <- pwc$min.p[pwc$stop.method == x$stopMethods[k]]	
		collapsedUnw <- pwc$min.p[pwc$stop.method == 'unw']		
		collRanks <- rank(collapsed, ties.method = "first")
		esBigHold <- collapsed > collapsedUnw
		esDat$tPVal[(1:(2*nVar)) + (cnt - 1)*2*nVar] <- c(collapsed, collapsedUnw)
		esDat$tRank[(1:(2*nVar)) + (cnt - 1)*2*nVar] <- rep(collRanks, 2)
	}
	}
   	
   	   #	if(is.null(subst))	subst <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	n.var2 <- max(esDat$tRank * (!is.na(esDat$tPVal)), na.rm=TRUE)
   	
   	pt1 <- xyplot(tPVal~tRank|whichComp, groups = weighted, data=esDat, xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "t- and chi-squared p-values \n (pairwise minimum)", pch = c(19,1), col = "black", scales = list(alternating = 1), par.settings = list(strip.background = list(col=stripBgCol)),
   	#subst = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subst]) & (esDat$tRank <= n.var2), 
   	ylim = c(-.1, 1.1), ..., 
   	   	panel = function(...){
   	   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col=ltBl, type="l")
#   		panel.abline(a=0, b=1, col="lightblue")
   		panel.xyplot(...)
   		}
)
   	
   	}
   
   if (plots =="ks" || plots ==5) {  ## ks plot
   	if(x$estimand =="ATE") pwc <- pairwiseComparison(x, collapse.to = "covariate")

   	esDat <- makePlotDat(x$psList[[1]], whichPlot = 5, subsetStopMeth = subset)
   	nVar <- length(makePlotDat(x$psList[[1]], whichPlot = 5, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE))
   	effSzList <- effSzListUnw <- matrix(NA, nrow = nVar, ncol = length(x$psList))   	
   	
   	cnt <- 0
   	for(k in subset){
   		cnt <- cnt + 1
   	for(j in 1:length(x$psList)){	
	x2 <- x$psList[[j]]
	effSzList[,j] <- makePlotDat(x2, whichPlot = 5, subsetStopMeth = k, yOnly = TRUE, incUnw = FALSE)
	effSzListUnw[,j] <- makePlotDat(x2, whichPlot = 5, subsetStopMeth = 1, yOnly = TRUE, incUnw = TRUE)[1:nVar]
	collapsed <- apply(effSzList, 1, max)
	if(x$estimand == "ATE") collapsed <- pwc$min.ks.pval[pwc$stop.method == x$stopMethods[k]]
	collapsedUnw <- apply(effSzListUnw, 1, max)
	if(x$estimand == "ATE") collapsedUnw <- pwc$min.ks.pval[pwc$stop.method == "unw"]	
	collRanks <- rank(collapsed, ties.method = "first")
	esBigHold <- collapsed > collapsedUnw
	esDat$ksPVal[(1:(2*nVar)) + (cnt-1)*2*nVar] <- c(collapsed, collapsedUnw)
	esDat$ksRank[(1:(2*nVar)) + (cnt-1)*2*nVar] <- rep(collRanks, 2)
	}
	}


   	
   	if(is.null(subst))	subst <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	n.var2 <- max(esDat$ksRank*(!is.na(esDat$ksPVal)), na.rm=TRUE)
   	
   	pt1 <- xyplot(ksPVal~ksRank|whichComp, groups=weighted, scales = list(alternating = 1), data = esDat,ylim = c(-.1, 1.1), ..., xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "KS test p-values", pch = c(19,1), col="black",par.settings = list(strip.background = list(col=stripBgCol)),
   	subst = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subst]) & (esDat$ksRank <= n.var2),
   	panel = function(...){
   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col=ltBl, type="l")
   		panel.xyplot(...)
   	})
   	}
   	
return(pt1)
	
	}

}

}
