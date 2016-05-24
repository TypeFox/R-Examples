##################################################################################
# plot function spotPlot() 
# collection of plotting procedures
###################################################################################

###################################################################################
#' Read .bst File
#'
#' This function reads the .bst file of the current project. The data is used
#' for plots and reports. Will usually not be used if no result/best files are written.
#' 
#' @param spotConfig the list of all parameters is given, used here 
#' 		\code{spotReadBstFile}: the name pf the .bst file to be read
#' 		\code{spotConfig$io.columnSep}: Separator for the columns 
#'
#' @return data.frame \code{bstData} \cr
#' - \code{bstData} holds a column "y" with the results and all columns with column-names 
#' derived from .roi file (should be the parameters of the algorithm
#'
#' @seealso \code{\link{SPOT}} \code{\link{spot}} \code{\link{spotPlotBst}}
#' @export
####################################################################################
spotReadBstFile<-function(spotConfig){
	read.table(spotConfig$io.bstFileName
			, sep=spotConfig$io.columnSep
			, header = TRUE	
			, stringsAsFactors = TRUE);
}

###################################################################################
#' Plot Best Solution Found so far
#'
#' Function used to continuously plot the actually retrieved best value 
#' throughout a SPOT run. The number of variables shown is limited to twelve -
#' make sure the relevant variables belong to the first lines of your .roi-file.
#' 
#' @param spotConfig the list of all parameters is given, used here 
#' 		\code{alg.roi}: the region of interest reduced to a matrix
#' @seealso \code{\link{SPOT}} \code{\link{spot}} \code{\link{spotReadBstFile}}
#' @export
#' @keywords internal
####################################################################################
spotPlotBst <- function(spotConfig){	
	rawB <- spotGetRawDataMatrixB(spotConfig);
	nEvals <- nrow(rawB)
	pNames <- row.names(spotConfig$alg.roi);
	if(spotConfig$spot.fileMode){
		b<-spotReadBstFile(spotConfig)	
	}else{
		b<-spotConfig$alg.currentBest
	}
	b<-unique(b)#fix to avoid multiple bests from 1 step. requirement is that steps are included in best file	
	Y = b[,spotConfig$alg.resultColumn]
	n <- max(spotConfig$alg.currentResult$STEP)+1
	step<-1:n
	# do not try to show more than 11 variables or the graphic window will crash #todo depends on machine
	numPlots<-min(length(pNames),11,spotConfig$report.io.maxPlots)
	if(numPlots==0) 
		return # if there is nothing to do, do not try to 
	if(numPlots>=4){ # more than 3 variables are not nice to plot, but possible: 
		nRows <- ceiling(sqrt(1+numPlots))
		nCols <- ceiling(numPlots/sqrt(1+numPlots))
 	} else {	## if we have just 3 variables, we can use a simple version here:
		nRows <- 1+numPlots
		nCols <- 1
	}	
	if(length(spotConfig$alg.resultColumn)==1){
		par(mfrow=c(nRows, nCols), xpd=NA)
		plot(step,Y, type="b", ylab = "Y", main = paste("Eval: ",as.character(nEvals), ", Y: ", as.character(Y[length(Y)])))
		for (i in 1:numPlots){
			plot(step,t(b[pNames[i]]), type="b", ylab = as.character(pNames[i]))	
		}
	}	
	else{
		nCols=2
		nRows=1
		par(mfrow=c(nRows, nCols), xpd=NA,oma = c(1.3, 0, 0, 0))
		mD<-spotPrepareData(spotConfig) 
		sel<-which(mD$STEP==max(mD$STEP))
		#first plot: last evaluated design
		plot(mD$mergedY[sel,1],mD$mergedY[sel,2], type="p",
			xlab=spotConfig$alg.resultColumn[1], ylab = spotConfig$alg.resultColumn[2],
			main = "Last sequential step")
		#second plot: dominating points in all evaluations
		sel<-which(!is_dominated(t(as.matrix(mD$mergedY))))
		plot(mD$mergedY[sel,1],mD$mergedY[sel,2], type="p",
			xlab=spotConfig$alg.resultColumn[1], ylab = spotConfig$alg.resultColumn[2],
			main = "Pareto front (all steps)")
		mtext(paste("Eval: ",as.character(nEvals),"   Step: ",as.character(max(mD$STEP))), side=1, outer = TRUE, cex = 1)			
	}
}


