diag.plot.color <- function(x,plots, subset, color, time = NULL, ...)
{
	ltBl <- ifelse(color, "lightblue","gray80")
	rdCol <- ifelse(color, "red","black")
	stripBgCol <- ifelse(color, "#ffe5cc", "transparent")
	ptSymCol <- ifelse(color, "#0080ff", "black")
	dots <- list(...)

	treat <- x$treat
	propScores <- x$ps
	weights <- x$w
   if (!all(weights[,1]==1)){
      weights   <- cbind(unw=rep(1,nrow(weights)),weights)
      propScores <- cbind(unw=rep(0.5,nrow(propScores)),propScores)
   }
    
whichVar <- pVal <- weighted <- varb <- NULL   

if (plots == "optimize" || plots == 1) {
	
	optDat <- makePlotDat(x, whichPlot = 1)

	if(is.null(subset))
	subset <- 1:length(levels(as.factor(optDat$stopRule)))
	
	if(is.null(time)){
		xlb <- "Iteration"
	}
	else{
		xlb <- paste("Iteration (Time ", time, ")", sep = "")
	}

	pt1 <- xyplot(balance ~ iteration | stopRule, data = optDat, ylab = "Balance measure", xlab = xlb, scales = list(alternating = 1), as.table = TRUE, subset = as.factor(optDat$stopRule) %in% levels(as.factor(optDat$stopRule))[subset], col = ptSymCol, par.settings = list(strip.background = list(col=stripBgCol)), ...)
	
}
   
   if (plots == "es" || plots == 3)	{ ## es plot
   	
   	esDat <- makePlotDat(x, whichPlot = 3, yOnly = FALSE)
   	
   	if(is.null(subset))
   	subset <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	yMax <- min(3,max(esDat$effectSize, na.rm=TRUE)) + .05	
   	
   	nullPlot <- TRUE
   	
   	if(max(esDat[,1], na.rm=TRUE) > 3)
   	warning("Some effect sizes are larger than 3 and may not have been plotted.\n")	
   	subsetHold <- !esDat$esBig & (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset])
   	
	if(is.null(time)){
		xlb <- NULL
	}
	else{
		xlb <- paste("Time ", time, sep = "")
	}   	
   	
   	if(any(subsetHold)){   	
   		esDatTmp <- esDat
   		esDatTmp$effectSize[!subsetHold] <- NA	
   		pt1.1 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, data = esDatTmp, scales = list(alternating = 1),
   		ylim = c(-.05, yMax), type = "l", col = ltBl, as.table = TRUE,
   		subset = subsetHold, par.settings = list(strip.background = list(col=stripBgCol)),
   		ylab = "Absolute standard difference", xlab = xlb, ...,
   		panel = function(...){
   			panel.abline(h=c(.2,.5,.8), col="gray80")
   			panel.xyplot(...)
   		
   		})
   		ptHold <- pt1.1
   		nullPlot <- FALSE
   	}
   	
   	subsetHold <- esDat$esBig & (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset])
   	#subsetHold <- as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset]
   	
   	if(any(subsetHold)){
   	   	esDatTmp <- esDat
   		esDatTmp$effectSize[!subsetHold] <- NA		
   		pt1.2 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, subset = subsetHold, 
   		data = esDatTmp, ylab = "Absolute standard difference", xlab = xlb, as.table = TRUE, 
   		ylim = c(-.05, yMax), type = "l", col = rdCol, par.settings = list(strip.background = list(col=stripBgCol)),
   		lwd = 2, ...)
   		if(nullPlot){
   			ptHold <- pt1.2
   			nullPlot <- FALSE
   		}
   		else {
   			ptHold <- ptHold + as.layer(pt1.2)
   			}
   	}
   	
   	subsetHold <- as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset]
   	
   	if(all(esDat$pVal < 0.05, na.rm=TRUE)) pchHold <- 19
   	else if(all(esDat$pVal >= 0.05, na.rm=TRUE)) pchHold <- 1
   	else pchHold <- c(19,1)
   	
   	pt2 <- xyplot(effectSize ~ weighted | whichComp, groups = pVal >= .05, data = esDat,
   		ylab = "Absolute standard difference", xlab = xlb, 
   		ylim = c(-.05, yMax), type = "p", col = rdCol, pch = pchHold,
   		subset = subsetHold, par.settings = list(strip.background = list(col=stripBgCol)), ...)
   	ptHold <- ptHold + pt2
   	
   	pt1 <- ptHold
   	
   	}
   						
 			
   
   if (plots == "t" || plots == 4) { ## t p-values plot
   	
   	if(is.null(time)){
		xlb <- "Rank of p-value for pretreatment variables \n (hollow is weighted, solid is unweighted)"
	}
	else{
		xlb <- paste("Rank of p-value for pretreatment variables \n (hollow is weighted, solid is unweighted) \n Time ", time, sep = "")
	}   	
   	
   	
   	esDat <- makePlotDat(x, whichPlot = 4)
   	if(is.null(subset)) subset <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	
   	
   	n.var2 <- max(esDat$tRank * (!is.na(esDat$tPVal)), na.rm=TRUE)
   	pt1 <- xyplot(tPVal~tRank|whichComp, groups = weighted, data=esDat, xlab = xlb, ylab = "T test p-values", pch = c(19,1), col = "black", scales = list(alternating = 1),
   	subset = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset]) & (esDat$tRank <= n.var2), ylim = c(-.1, 1.1), par.settings = list(strip.background = list(col=stripBgCol)),..., 
   	   	panel = function(...){
   	   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col=ltBl, type="l")
	   		panel.xyplot(...)
   		}
)
   	
   	}
   
   if (plots =="ks" || plots ==5) {  ## ks plot
   	esDat <- makePlotDat(x, whichPlot = 5)
   	
   	if(is.null(subset)) subset <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	n.var2 <- max(esDat$ksRank*(!is.na(esDat$ksPVal)), na.rm=TRUE)
   	
   	if(is.null(time)){
		xlb <- "Rank of p-value for pretreatment variables \n (hollow is weighted, solid is unweighted)"
	}
	else{
		xlb <- paste("Rank of p-value for pretreatment variables \n (hollow is weighted, solid is unweighted) \n Time ", time, sep = "")
	}   	   	
   	
   	pt1 <- xyplot(ksPVal~ksRank|whichComp, groups=weighted, scales = list(alternating = 1), data = esDat,ylim = c(-.1, 1.1), ..., xlab = xlb, ylab = "KS p-values", pch = c(19,1), col="black", par.settings = list(strip.background = list(col=stripBgCol)),
   	subset = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset]) & (esDat$ksRank <= n.var2),
   	panel = function(...){
   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col= ltBl, type="l")
   		panel.xyplot(...)
   	})
   	}
   	
#   	if (plots == "histogram" || plots == 10)
#   	pt1 <- histogram.dxwts(...)
   	
   	if (plots == "boxplot" || plots == 2){
			pt1 <- boxplot(x, color = color, subset = subset, time = time, ...)
   		}
   	
   	if (plots == "histogram" || plots == 6){
   		
   		treat <- x$treat
		wghts <- x$w
		vars <- x$variableNames
		if (is.null(dots$main)) dots$main <- "Control weights"
#	if(is.null(x$w.ctrl) || is.null(x$treat))
#         stop("For the weight histogram w.ctrl and treat cannot be NULL")	
	if(!all(wghts[,1]==1))   {
		controlWeights <- wghts[x$treat == 0, ]
	}
	else controlWeights <- wghts[x$treat == 0, -1]
	
	controlWeights <- as.matrix(controlWeights)
	nFrame <- ncol(controlWeights)
	longWeights <- matrix(t(controlWeights), ncol=1)
	longWeights <- data.frame(Weights = longWeights, varb = names(wghts))
	
	if(is.null(subset))
	subset <- 1:length(unique(longWeights$varb))
	
	if(color == FALSE & is.null(dots$col)) 	
	pt1 <- histogram(~Weights|varb, data=longWeights, par.settings = list(strip.background = list(col=stripBgCol)), subset = varb %in% unique(controlWeights)[subset], col = NULL, ...)
	else pt1 <- histogram(~Weights|varb, data=longWeights, par.settings = list(strip.background = list(col=stripBgCol)), subset = varb %in% unique(longWeights$varb)[subset], ...) 

   		
   		
   	}
   	
   	
   	if(!(plots %in% c(1:6, "boxplot","t","ks","optimize","es", "histogram")))
   	stop("plots must be an integer from 1 to 6, or be one of 'optimize' \n
   	'boxplot','es','t', 'ks', or 'histogram'\n")
   	
   
   return(pt1)


}
