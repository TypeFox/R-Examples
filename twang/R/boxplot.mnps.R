boxplot.mnps <- function(x, stop.method = NULL, color = TRUE, figureRows = NULL, singlePlot = NULL, multiPage = FALSE, time = NULL, print = TRUE,  ...){
	
	ptSymCol <- ifelse(color, "#0080ff", "black")	
	bwCols <- list(col = ptSymCol)
	stripBgCol <- ifelse(color, "#ffe5cc", "transparent")
	nPlot <- x$nFits
	
	
	ptHld <- vector(mode = "list", length = nPlot)
	
	
	if(is.null(stop.method)) stop.method <- x$stopMethods
	
	if(!is.null(singlePlot)){
		if(!is.numeric(singlePlot) | singlePlot > x$nFits) stop(paste("If specified, the \"singlePlot\" argument must be an integer between 1 and", x$nFits, "for this object."))
		if(round(singlePlot) != singlePlot) stop("If specified, the \"singlePlot\" argument must be a positive integer.")
	}
	
	if(length(stop.method) > 1){ 
		if(is.numeric(stop.method))
		warning("Using only the first stop.method, ",  x$stopMethods[1])
		else warning("Using only the first stop.method, ",  stop.method[1])
		stop.method <- stop.method[1]
		}
			
	
	
	if(is.numeric(stop.method)) stop.method = x$stopMethods[stop.method]
	

	
	if(x$estimand == "ATE"){

		bwDat <- whichResp <- NULL
			
		stopMethLong <- paste(stop.method, ".ATE", sep = "")
			
			
		for(j in 1:nPlot){
			bwDat <- data.frame(ps = x$psList[[j]]$ps[,stopMethLong], treat = x$treatVar, whichResp = x$treatLev[j])
			ptNm <- paste(x$treatLev[j], " propensity scores by Tx group")
			if(!is.null(time)) ptNm <- paste(ptNm, " (time ", time, ")", sep = "")
			pt1 <- bwplot(ps ~ treat, groups = whichResp,  
		xlab = "Treatment", ylab = "Propensity scores", ylim = c(-.1,1.1), data = bwDat, main = ptNm, par.settings = list(strip.background = list(col=stripBgCol), box.rectangle = bwCols, plot.symbol = bwCols, box.umbrella = bwCols), ...)
		
		ptHld[[j]] <- pt1

		}
		}
		
		
	else if(x$estimand == "ATT"){
		
		bwDat <- NULL
		
		stopMethLong <- paste(stop.method, ".ATT", sep = "")
		bwDat <- NULL
		
		for(j in 1:nPlot){
			currCats <- c(x$treatATT, x$levExceptTreatATT[j])
			bwDat <- data.frame(ps = x$psList[[j]]$ps[,stopMethLong], treat = currCats[1 + x$psList[[j]]$data$currResp], respCat = x$levExceptTreatATT[j], attGrp = x$treatATT)
			ptNm <- paste("Propensity score of ", x$levExceptTreatATT[j], " versus ", x$treatATT, ".", sep = "")
			if(!is.null(time)) ptNm <- paste(ptNm, " (time ", time, ")", sep = "")
			pt1 <- bwplot(ps ~ treat, data = bwDat, ylim = c(-.1,1.1), ylab = "Propensity scores", xlab = "Treatment", main = ptNm,par.settings = list(strip.background = list(col=stripBgCol), box.rectangle = bwCols, plot.symbol = bwCols, box.umbrella = bwCols), ...)
			
			ptHld[[j]] <- pt1
	
			
		}		
	}	
	
	if(print) displayPlots(ptHld, figureRows = figureRows, singlePlot = singlePlot, multiPage = multiPage, bxpt = TRUE)
	return(ptHld)


}