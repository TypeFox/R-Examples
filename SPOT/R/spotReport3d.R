
# TODO: force model and observations to same scaling, if desired

###################################################################################################
#' 3d Plot of Meta Model - Report Function
#' 
#' Function to generate a 3d surface plot of the predicted meta model.
#'
#' This report function uses the parameter spotConfig$seq.modelFit to plot the predicted model.
#' If spotConfig$seq.modelFit is NULL, the model is generated, based on spotConfig$seq.predictionModel.func.
#' It is not recommended to use this function at the end of an "auto" run of SPOT, make sure to save results first. 
#' By default, twiddler will be used to let the user specify which of the parameters should be varied in the plot.
#' Values not varied in the graph are fixed to their "best" value according to the current Results.
#'
#' @param spotConfig the configuration list of all spot parameters. \cr
#'				The parameter spotConfig$report.interactive=TRUE will be set as default if not contained in the list.
#'				That means, by default the user will be asked to specify which parameters will be varied when the report is started. This is done in a small twiddler GUI. 
#'				If the user wants to specify which parameters should be plotted against each other, before starting the report, he can 
#'				set the parameter spotConfig$report.aIndex and spotConfig$report.bIndex. They should be two different integer numbers. They will only be used if 
#'				spotConfig$report.interactive is FALSE. By default they will be set to 1 and 2, so the first two parameters in the ROI will be plotted.
#'				In case of a multi objective problem, spotConfig$report.cIndex will determine which point from the Pareto front will be used to determine values of 
#'				parameters not on the axis of the plot.
#' @seealso \code{\link{spotSurfContour}}
#' @export
###################################################################################################
spotReport3d <- function(spotConfig) {
	#set unset defaults
	if(is.null(spotConfig$report.interactive)){spotConfig$report.interactive=TRUE}
	if(is.null(spotConfig$report.main)){spotConfig$report.main="predicted Model"}
	if(is.null(spotConfig$report.3d.clipZ)){spotConfig$report.3d.clipZ=NA}
	if(is.null(spotConfig$report.observations)){spotConfig$report.observations=TRUE}
	#in case of MCO
	if(length(spotConfig$alg.resultColumn)>1){ #generate a plot for each objective ....
		spotConfig<-spotReport3dMulti(spotConfig)
		return(spotConfig)
	}
	#load packages
	spotInstAndLoadPackages("rgl")
	spotInstAndLoadPackages("twiddler")		
	spotWriteLines(spotConfig$io.verbosity,2,"  Creating 3d-report:");	
	rawB <- spotGetRawDataMatrixB(spotConfig);
	spotPrint(spotConfig$io.verbosity,1,summary(rawB));
	mergedData <- spotPrepareData(spotConfig)
	mergedB <- spotGetMergedDataMatrixB(mergedData, spotConfig);	
	xNames <- row.names(spotConfig$alg.roi)
	yNames <- spotConfig$alg.resultColumn
	spotConfig=spotWriteBest(mergedData, spotConfig);
	C1=spotConfig$alg.currentBest[nrow(spotConfig$alg.currentBest),]
	spotWriteLines(spotConfig$io.verbosity,1," ");
	spotPrint(spotConfig$io.verbosity,1,paste("Best solution found with ",nrow(rawB)," evaluations:",sep=""));
	spotPrint(spotConfig$io.verbosity,1,C1);	
	a=1;
	b=2;
	plotFn <- function(aIndex, bIndex) {
		if((aIndex==bIndex)){#here: add more errors
			warning("Plot index problem: both chosen indices are the same") 
		}
		lo <- spotConfig$alg.roi$lower[c(aIndex,bIndex)]
		up <- spotConfig$alg.roi$upper[c(aIndex,bIndex)]
		fn <- function(a,b){
			x<-matrix(t(C1[xNames]),nrow(spotConfig$alg.roi),length(a))
			x<-t(x)
			x[,aIndex]=a
			x[,bIndex]=b
			colnames(x) <- row.names(spotConfig$alg.aroi)
			if(is.null(spotConfig$seq.modelFit)){
				spotConfig1 <- eval(call(spotConfig$seq.predictionModel.func
                                        , rawB
                                        , mergedB
                                        , as.data.frame(x)
                                         ,spotConfig))
			}else{
				spotConfig1 <- eval(call(spotConfig$seq.predictionModel.func
                                , NULL 
								, NULL
								, as.data.frame(x)
								, spotConfig
                                , spotConfig$seq.modelFit #external fit is used, model is only evaluated not build
								))
			}
			spotConfig1$seq.largeDesignY[[1]];
		}
		clip = spotConfig$report.3d.clipZ
		if(!any(is.na(clip)) && length(clip)==1){
			clip = range(rawB[yNames])
		}
		if(spotConfig$report.observations){
			datx<-mergedB[,rownames(spotConfig$alg.roi)]
			dat<-data.frame(datx[,c(aIndex,bIndex)],mergedB[yNames])
			spotSurf3d(fn,lo,up,100,clip,dat)
		}else{spotSurf3d(fn,lo,up,100,clip)}
	}
	if(spotConfig$report.interactive==TRUE){
		twiddle(plotFn(a, b), eval = FALSE, a = knob(c(1, nrow(spotConfig$alg.roi)), res = 1),
			b = knob(c(1, nrow(spotConfig$alg.roi)), res = 1, default=2))
	}
	else{
		if(is.null(spotConfig$report.aIndex)){spotConfig$report.aIndex=1}
		if(is.null(spotConfig$report.bIndex)){spotConfig$report.bIndex=2}
		plotFn(spotConfig$report.aIndex,spotConfig$report.bIndex)
	}	
	spotConfig
}

###################################################################################################
#' 3d Plot of Meta Model - Report Function
#' 
#' Function to generate a 3d surface plot of the predicted meta model.
#'
#' This report function uses the parameter spotConfig$seq.modelFit to plot the predicted model.
#' If spotConfig$seq.modelFit is NULL, the model is generated, based on spotConfig$seq.predictionModel.func.
#' It is not recommended to use this function at the end of an "auto" run of SPOT, make sure to save results first. 
#' By default, twiddler will be used to let the user specify which of the parameters should be varied in the plot.
#' Values not varied in the graph are fixed to their "best" value according to the current Results.
#'
#' @param spotConfig the configuration list of all spot parameters. \cr
#'				The parameter spotConfig$report.interactive=TRUE will be set as default if not contained in the list.
#'				That means, by default the user will be asked to specify which parameters will be varied when the report is started. This is done in a small twiddler GUI. 
#'				If the user wants to specify which parameters should be plotted against each other, before starting the report, he can 
#'				set the parameter spotConfig$report.aIndex and spotConfig$report.bIndex. They should be two different integer numbers. They will only be used if 
#'				spotConfig$report.interactive is FALSE. By default they will be set to 1 and 2, so the first two parameters in the ROI will be plotted.
#' @export
#' @keywords internal
###################################################################################################
spotReport3dMulti <- function(spotConfig) {
	#load packages
	spotInstAndLoadPackages("rgl")
	spotInstAndLoadPackages("twiddler")		
	spotWriteLines(spotConfig$io.verbosity,2,"  Creating 3d-report for multiple objectives:")
	rawB <- spotGetRawDataMatrixB(spotConfig)
	spotPrint(spotConfig$io.verbosity,1,summary(rawB))
	mergedData <- spotPrepareData(spotConfig)
	mergedB <- spotGetMergedDataMatrixB(mergedData, spotConfig);	
	yNames <- spotConfig$alg.resultColumn
	## calculate pareto front and set and hvolume
	allPoints<-t(mergedData$mergedY)
	frontPoints<-nondominated_points(allPoints)
	if(is.null(spotConfig$mco.refPoint)){
		vol<-dominated_hypervolume(frontPoints)
	}else{
		vol<-dominated_hypervolume(frontPoints,ref=spotConfig$mco.refPoint)
	}	
	spotPrint(spotConfig$io.verbosity,1,paste("Hypervolume indicator of pareto optimal solutions: ",vol,sep=""))
	spotConfig$mco.hvolume = vol # hyper volume
	spotConfig$mco.val= t(frontPoints) #pareto front
	spotConfig$mco.par=  mergedData$x[!is_dominated(allPoints),] #pareto set
	##
	spotConfig=spotWriteBest(mergedData, spotConfig)
	C1=spotConfig$alg.currentBest[nrow(spotConfig$alg.currentBest),]
	spotWriteLines(spotConfig$io.verbosity,1," ")
	spotPrint(spotConfig$io.verbosity,1,paste("Best solution found with ",nrow(rawB)," evaluations:",sep=""))
	spotPrint(spotConfig$io.verbosity,1,C1)
	a=1;
	b=2;
	plotFn <- function(aIndex, bIndex, cIndex) {
		if((aIndex==bIndex)){#here: add more errors
			warning("Plot index problem: both chosen indices are the same") 
		}
		lo <- spotConfig$alg.roi$lower[c(aIndex,bIndex)]
		up <- spotConfig$alg.roi$upper[c(aIndex,bIndex)]
		fn <- function(a,b,ii){
			x<-matrix(t(spotConfig$mco.par[cIndex,]),nrow(spotConfig$alg.roi),length(a))
			x<-t(x)
			x[,aIndex]=a
			x[,bIndex]=b
			colnames(x) <- row.names(spotConfig$alg.roi);	
			if(is.null(spotConfig$seq.modelFit)){
				spotConfig1 <- eval(call(spotConfig$seq.predictionModel.func
                                        , rawB
                                        , mergedB
                                        , as.data.frame(x)
                                         ,spotConfig));
			}else{
				spotConfig1 <- eval(call(spotConfig$seq.predictionModel.func
                                , NULL 
								, NULL
								, as.data.frame(x)
								, spotConfig
                                , spotConfig$seq.modelFit #external fit is used, model is only evaluated not build
								));
			}
			spotConfig1$seq.largeDesignY[[ii]]
		}
		for(i in 1:length(spotConfig$alg.resultColumn)){
			clip = spotConfig$report.3d.clipZ
			if(!any(is.na(clip)) && length(clip)==1){
				clip =  range(rawB[yNames[i]])
			}
			if(spotConfig$report.observations){
				datx<-mergedB[,rownames(spotConfig$alg.roi)]
				dat<-data.frame(datx[,c(aIndex,bIndex)],mergedB[yNames][,i])
				spotSurf3d(fn,lo,up,100,clip,dat,ii=i)
			}else{spotSurf3d(fn,lo,up,100,clip,ii=i)}
			#			
		}	
	}
	if(spotConfig$report.interactive==TRUE){
		twiddle(plotFn(a, b, C), eval = FALSE, a = knob(c(1, nrow(spotConfig$alg.roi)), res = 1),
			b = knob(c(1, nrow(spotConfig$alg.roi)), res = 1, default=2),
			C = knob(c(1, nrow(spotConfig$mco.par)), res = 1))
	}
	else{
		if(is.null(spotConfig$report.aIndex)){spotConfig$report.aIndex=1}
		if(is.null(spotConfig$report.bIndex)){spotConfig$report.bIndex=2}
		if(is.null(spotConfig$report.bIndex)){spotConfig$report.cIndex=1}	
		plotFn(spotConfig$report.aIndex,spotConfig$report.bIndex,spotConfig$report.cIndex)
	}	
	spotConfig
}