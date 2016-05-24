diag.plot <- function(x,plots, subset, ...)
{
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

	pt1 <- xyplot(balance ~ iteration | stopRule, data = optDat, ylab = "Balance measure", xlab = "Iteration", scales = list(alternating = 1), subset = as.factor(optDat$stopRule) %in% levels(as.factor(optDat$stopRule))[subset], ...)	
	
	
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
   	
   	if(any(subsetHold)){
   	pt1.1 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, data = esDat, scales = list(alternating = 1),
   	ylim = c(-.05, yMax), type = "l", col = "lightblue", 
   	subset = subsetHold, 
   	ylab = "Absolute standard difference", xlab = NULL, ...,
   	panel = function(...){
   		panel.abline(h=c(.2,.5,.8), col="gray80")
   		panel.xyplot(...)
   		
   	})
   	ptHold <- pt1.1
   	nullPlot <- FALSE
   	}
   	
   	subsetHold <- esDat$esBig & (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset])
   	
   	if(any(subsetHold)){
   	pt1.2 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, 
   	data = esDat, ylab = "Absolute standard difference", xlab = NULL,
   	ylim = c(-.05, yMax), type = "l", col = "red", 
   	subset = subsetHold, lwd = 2)
   	if(nullPlot){
   		ptHold <- pt1.2
   		nullPlot <- FALSE
   	}
   	else ptHold <- ptHold + pt1.2
   	}
   	
   	subsetHold <- as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset]
   	
   	if(all(esDat$pVal < 0.05, na.rm=TRUE)) pchHold <- 19
   	else if(all(esDat$pVal >= 0.05, na.rm=TRUE)) pchHold <- 1
   	else pchHold <- c(19,1)
   	
   	if(any(subsetHold)){
   	pt2 <- xyplot(effectSize ~ weighted | whichComp, groups = pVal, data = esDat,
   	ylab = "Absolute standard difference", xlab = NULL, 
   	ylim = c(-.05, yMax), type = "p", col = "red", pch = pchHold,
   	subset = subsetHold)
   	ptHold <- ptHold + pt2
   	}
   	
   	pt1 <- ptHold
   	
   	}
   						
 			
   
   if (plots == "t" || plots == 4) { ## t p-values plot
   	
   	esDat <- makePlotDat(x, whichPlot = 4)
   	if(is.null(subset)) subset <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	
   	
   	n.var2 <- max(esDat$tRank * (!is.na(esDat$tPVal)), na.rm=TRUE)
   	
   	pt1 <- xyplot(tPVal~tRank|whichComp, groups = weighted, data=esDat, xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "T test p-values", pch = c(19,1), col = "black", scales = list(alternating = 1),
   	subset = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset]) & (esDat$tRank <= n.var2), ylim = c(-.1, 1.1), ..., 
   	   	panel = function(...){
   	   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col="lightblue", type="l")
#   		panel.abline(a=0, b=1, col="lightblue")
   		panel.xyplot(...)
   		}
)
   	
   	}
   
   if (plots =="ks" || plots ==5) {  ## ks plot
   	esDat <- makePlotDat(x, whichPlot = 5)
   	
   	if(is.null(subset))
   	subset <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	n.var2 <- max(esDat$ksRank*(!is.na(esDat$ksPVal)), na.rm=TRUE)
   	
   	pt1 <- xyplot(ksPVal~ksRank|whichComp, groups=weighted, scales = list(alternating = 1), data = esDat,ylim = c(-.1, 1.1), ..., xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "KS test p-values", pch = c(19,1), col="black",
   	subset = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset]) & (esDat$ksRank <= n.var2),
   	panel = function(...){
   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col="lightblue", type="l")
   		panel.xyplot(...)
   	})
   	}
   	
#   	if (plots == "histogram" || plots == 10)
#   	pt1 <- histogram.dxwts(...)
   	
   	if (plots == "boxplot" || plots == 2)
   	pt1 <- boxplot(x, ...)
   	
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
	longWeights <- data.frame(Weights = longWeights, varb = colnames(controlWeights))
	
	if(is.null(subset))
	subset <- 1:length(colnames(controlWeights))
	
	
	pt1 <- histogram(~Weights|varb, data=longWeights, subset = varb %in% colnames(controlWeights)[subset],...)

   		
   		
   	}
   	
   	
   	if(!(plots %in% c(1:6, "boxplot","t","ks","optimize","es", "histogram")))
   	stop("plots must be an integer from 1 to 5, or be one of 'optimize' \n
   	'boxplot','es','t', 'ks', or 'histogram'\n")
   	
   
   return(pt1)


}
