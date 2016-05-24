plot.iptw <- function(x,plots="optimize", subset = NULL, color = TRUE, timePeriods = NULL, multiPage = FALSE, figureRows = NULL, ...)
{
   # Creates diag.plot plots and sends to current device
   # x:     ps object 
   # label: Label added to the plot titles

   # extract the propensity scores and weights from the ps object
   
   singlePlot <- NULL
   
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
   
   if(is.null(timePeriods)) timePeriods <- 1:x$nFits
   
   hdPt <- vector(mode = "list", length = length(timePeriods))
   
   for(i in 1:length(timePeriods)){
   	hdPt[[i]] <- diag.plot.color(x$psList[[i]], plots, subset = subset, color = color, time = timePeriods[i], ...)
   }
   
   	if(length(timePeriods) == 1) return(hdPt[[1]])
	else displayPlots(hdPt, figureRows = figureRows, singlePlot = singlePlot, multiPage = multiPage)	
}
