plotTruePairs <- function(x, subset = NULL, color = TRUE, ...){

   	ltBl <- ifelse(color, "lightblue","gray80")
	rdCol <- ifelse(color, "red","black")
	stripBgCol <- ifelse(color, "#ffe5cc", "transparent")
	ptSymCol <- ifelse(color, "#0080ff", "black")
	
	stop.method <- whichVar <- sig <- NULL
	
	
	pairs <- pairwiseComparison(x)
	stp <- x$stopMethods
	unwTab <- subset(pairs, stop.method == "unw")
	hdUnwTab <- unwTab
	if(length(stp) > 1){
		for(i in 1:(length(stp) - 1)) unwTab <- rbind(unwTab, hdUnwTab)
	}
	
	wghtTab <- subset(pairs, stop.method != "unw")
	
	unwTab <- data.frame(unwTab, whichVar = 1:nrow(unwTab))
	wghtTab <- data.frame(wghtTab, whichVar = 1:nrow(wghtTab))
	
	esDat <- rbind(unwTab, wghtTab)
	
	esDat$stopMeth <- wghtTab$stop.method
	
	esDat$weighted <- factor(rep(c("Unweighted","Weighted"), each = (nrow(esDat)/2)))
	
	esDat$bigger <- rep(as.numeric(wghtTab$std.eff.sz > unwTab$std.eff.sz),2)
	
	esDat$sig <- as.numeric(esDat$p < 0.05)
	
	esDat <- esDat[,c("tmt1","tmt2","whichVar","std.eff.sz","p","stopMeth","bigger","sig", "weighted")]
	
	names(esDat) <- c("tmt1","tmt2","whichVar","effectSize","pVal","whichComp","esBig","sig", "weighted")
	
	esDat$effectSize <- abs(esDat$effectSize)

   	yMax <- min(3,max(esDat[,"effectSize"])) + .05	
   	
   	if(max(esDat[,"effectSize"], na.rm=TRUE) > 3)
   	warning("Some effect sizes are larger than 3 and may not have been plotted.\n")	
   	
    
    nullPlot <- TRUE
       	
   	subsetHold <- !esDat$esBig 
   	
   	
   	if(any(subsetHold)){
   		esDatTmp <- esDat
   		esDatTmp$effectSize[!subsetHold] <- NA
   		pt1.1 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, data = esDatTmp, scales = list(alternating = 1),
   			ylim = c(-.05, yMax), type = "l", col = ltBl, as.table = TRUE,
   			ylab = "Absolute standardized difference \n (all pairs)", xlab = NULL, par.settings = list(strip.background = list(col=stripBgCol)),
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
   	
   	
   	if(all(esDat$sig >= .05)) pchHold <- 19
   	else if(all(esDat$sig < .05)) pchHold <- 1
   	else pchHold <- c(19,1)
   	
   	pt2 <- xyplot(effectSize ~ weighted | whichComp, groups = (sig < 0.05), data = esDat,
   	ylab = "Absolute standard difference", xlab = NULL, as.table = TRUE, 
   	ylim = c(-.05, yMax), type = "p", col = rdCol, pch = pchHold,par.settings = list(strip.background = list(col=stripBgCol)),
   	...)
   	if(!nullPlot) currPt <- currPt + pt2
   	else currPt <- pt2
   	
   	return(currPt)
}