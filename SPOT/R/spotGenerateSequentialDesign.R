############################################################################################
#' Generate Design for next sequential evaluation
#' 
#' Creates a new design. Design points are determined with respect to the current result file. 
#' 
#' Uses the functions \link{spotPrepareData},\link{spotGetRawDataMatrixB},\link{spotGetMergedDataMatrixB},\link{spotWriteLines}
#' \link{spotWriteBest},\link{spotPlotBst}
#' returns a sequential design to be written to the file <xxx>.des (will be determined from calling function)
#' 
#' @param spotConfig the list of all parameters is given, but the used ones are: \cr
#'   spotConfig$io.resFileName: is checked for existence. If not found, function fails with error\cr
#'   spotConfig$algSourceSrcPath: needed for the error message \cr
#'   spotConfig$userConfFileName: needed for the error message \cr
#'   spotConfig$io.verbosity: needed for command window output \cr
#' 
#' @return data.frame \code{design} \cr
#' - \code{design} contains one or more new design points to be calculated 
#' @export
#' @keywords internal
############################################################################################
spotGenerateSequentialDesign <- function(spotConfig) {
	spotWriteLines(spotConfig$io.verbosity,2,"Entering generateSequentialDesign");	
	rawB <- spotGetRawDataMatrixB(spotConfig);
	mergedData <- spotPrepareData(spotConfig)
	mergedB <- spotGetMergedDataMatrixB(mergedData, spotConfig);
	
	#for continuing runs without saving results in spotConfig:
	if(is.null(spotConfig$alg.currentResult))spotConfig$alg.currentResult<- spotGetRawResData(spotConfig)$rawD;

	spotConfig=spotWriteBest(mergedData, spotConfig);
	if(spotConfig$io.verbosity>2){
		spotPlotBst(spotConfig)
	}
	
	#####################################################
	#####################################################
	#####################################################
	## (0) First check if ocba or linear budget allocation should be done
	## check for var of the y-values. Only if y varies (i.e. function is 
	### noisy and evaluations are repeated) use ocba. else create error msg.
	if(spotConfig$spot.ocba == TRUE){ # only needs to be checked for ocba=TRUE
		varies=TRUE;
		varY <- mergedData$varY
		for (i in 1:length(varY)){
			if (is.na(varY[i])||is.nan(varY[i])||is.null(varY[i])||(varY[i]==0)){
				varies=FALSE;
			}
		}
		if(varies==FALSE){
			stop("			
There is no variance for some point(s) in the current design.
Therefore OCBA cannot be used. Possible reasons are a target 
function without noise, or the design points are not repeated.
SPOT with OCBA makes only sense if the target function is noisy.
If a non noisy function is used, the default settings should 
be adopted,	as described in the help of spot() or spotOptim().
That means: either use spot.ocba=FALSE, or set the repeats
(init.design.repeats) to values larger
than 1.

The current variance vector for the design points is: 
",paste(varY," "))
		}
	}	
	

	#####################################################
	#####################################################
	#####################################################
	## (1) Here it is most important to cover a broad area of 
	## the search space, so Latin hypercube designs are preferred to factorial designs.
	## The user can specify what ever he wants...
	if(!exists(spotConfig$seq.design.func))stop(paste("The design function name", spotConfig$seq.design.func, "is not found in the workspace \n
		Please make sure to load the design function in the workspace, or specify the correct function in spotConfig$seq.design.func"))
	if(spotConfig$seq.design.size==0)
		largeDesign <- NULL
	else
		largeDesign <- (eval(call(spotConfig$seq.design.func, 
								spotConfig, 
								spotConfig$seq.design.size, 
								spotConfig$seq.design.retries)));
	#print(largeDesign)
			
	#####################################################
	#####################################################
	#####################################################
	### (2) Fit the prediction model and generate new sample points:
	### x contains input, y output values
	### now calling the seq.predictionModel.func specified in spotConfigure
	### the prediction model is build with the values from the resfiles
	if(!exists(spotConfig$seq.predictionModel.func))stop(paste("The prediction model function name", spotConfig$seq.predictionModel.func, "is not found in the workspace \n
		Please make sure to load the prediction model function in the workspace, or specify the correct function in spotConfig$seq.predictionModel.func" ))		
	spotConfig <- eval(call(spotConfig$seq.predictionModel.func
                                        , rawB
                                        , mergedB
                                        , largeDesign
                                         , spotConfig));

	#####################################################
	#####################################################	
	#####################################################
    ## (2b) If desired, optimize fit(s) returned by prediction model(s)
	if (!is.na(spotConfig$seq.predictionOpt.func)){
		spotConfig <- eval(call(spotConfig$seq.predictionOpt.func
											, largeDesign #for start point of optimization	
											, spotConfig));
		
		if(is.null(largeDesign)){
			largeDesign <- as.data.frame(spotConfig$optDesign)
			spotConfig$seq.largeDesignY	<- as.data.frame(spotConfig$optDesignY)	
		}else{
			largeDesign <-  as.data.frame(rbind(spotConfig$optDesign, as.matrix(largeDesign)),row.names=NULL); 
			spotConfig$seq.largeDesignY <-  as.data.frame(rbind(spotConfig$optDesignY, as.matrix(spotConfig$seq.largeDesignY)),row.names=NULL);  #TODO: check where largeDesignY is used, and if it is used correctly!
		}
		spotConfig$optDesignY<-NULL
		spotConfig$optDesign<-NULL
	}
	names(largeDesign)<- row.names(spotConfig$alg.roi);	
	#names(spotConfig$seq.largeDesignY)<-spotConfig$alg.resultColumn

	mcosort=TRUE
	if(is.null(ncol(spotConfig$seq.largeDesignY))){
		mcosort=FALSE
	}else{
		mcosort=if(ncol(spotConfig$seq.largeDesignY)==1){FALSE}else{TRUE}	
	}	

	#now sort the largeDesign, to select only best points for next evaluation
	if(length(spotConfig$alg.resultColumn)>1 & mcosort){ # only do this if mco (first condition) and model returns multiple objectives (second). Second will not be TRUE if a multi criteria based expected improvement is used, because this merges objectives into one
	#	in case of multi criteria spot: "sort" large design by hypervolume contribution nds rank
		#if(spotConfig$seq.mco.infill=="sort"){
		#	largeDesign <- spotMcoSort(largeDesign,spotConfig$seq.largeDesignY,spotConfig$seq.design.new.size)} #TODO only largedesign not largedesignY is sorted. For future applications might be needed to sort both
		#if(spotConfig$seq.mco.infill=="fill"){
			mergY <-  eval(call(spotConfig$seq.predictionModel.func, NULL, NULL  #reevaluate known points on model, to be in the same scale as the largeDesignY
								, mergedData$x
								, spotConfig
                                , spotConfig$seq.modelFit 
								))$seq.largeDesignY
								
			if(spotConfig$seq.mco.selection=="hypervol"){					
				largeDesign <- spotMcoSelectionHypervol(largeDesign,spotConfig$seq.largeDesignY,spotConfig$seq.design.new.size, mergedData$x,mergY,spotConfig$mco.refPoint)}
			if(spotConfig$seq.mco.selection=="tournament1"){	
				warning("Tournament selection with tournament1 is not recommended. Use hypervol or tournament2 for better results.")
				if(is.null(spotConfig$seq.mco.tournsize)){spotConfig$seq.mco.tournsize=1}				
				largeDesign <- spotMcoCrowdTournament(largeDesign,spotConfig$seq.largeDesignY,spotConfig$seq.mco.tournsize,spotConfig$seq.design.new.size)
			}
			if(spotConfig$seq.mco.selection=="tournament2"){					
				if(is.null(spotConfig$seq.mco.tournsize)){spotConfig$seq.mco.tournsize=1}				
				largeDesign <- spotMcoTournament(largeDesign,spotConfig$seq.largeDesignY,spotConfig$seq.mco.tournsize,spotConfig$seq.design.new.size)
			}		
		#}
	}else{ #in case of single criteria spot: sort large design by criteria value
		largeDesign <-  as.data.frame(largeDesign[order(spotConfig$seq.largeDesignY,decreasing=FALSE),]);
		#spotConfig$seq.largeDesignY <-  as.data.frame(spotConfig$seq.largeDesignY[order(spotConfig$seq.largeDesignY,decreasing=FALSE),]);
	}	
	maxpoints <- min(nrow(largeDesign),spotConfig$seq.design.new.size)
	largeDesignEvaluated <- as.data.frame(largeDesign[1:maxpoints,,drop=FALSE]); #limit to set design size, or to maximum available

    spotPrint(spotConfig$io.verbosity,3,"largeDesignEvaluated:")
	spotPrint(spotConfig$io.verbosity,3,largeDesignEvaluated)

	#####################################################
	#####################################################
	#####################################################
    ## (3) Adaptation of the number of repeats and
    ## (4) Combination of old (which should be re-evaluated)  and new design points 
	if (spotConfig$spot.ocba == TRUE){
		design <- spotRepeatsOcba(spotConfig,mergedData,largeDesignEvaluated)
	}else{		
		design <- spotRepeats(spotConfig,mergedData,largeDesignEvaluated)
	}	
	### finally write some screen output
	spotPrint(spotConfig$io.verbosity,2,"design:")
	spotPrint(spotConfig$io.verbosity,2,design)	
	spotWriteLines(spotConfig$io.verbosity,2,"  Leaving generateSequentialDesign");
	## write the design to the .des-file	
	if (spotConfig$spot.fileMode){
		spotWriteDes(design,spotConfig$io.verbosity,spotConfig$io.columnSep,spotConfig$io.desFileName)	
	}
	spotConfig$alg.currentDesign<-design;		
	spotConfig
}
