##################################################################################
#' Select repeats of design points based on OCBA
#'
#' @param spotConfig List, containing all settings of SPOT
#' @param mergedData merged data from runs of the target function
#' @param largeDesignEvaluated sorted large design in sequential SPOT step
#' 
#' @seealso  \code{\link{spotGenerateSequentialDesign}} \code{\link{spotOcba}} \code{\link{spotRepeats}}
#' @keywords internal
###################################################################################
spotRepeatsOcba <- function(spotConfig,mergedData,largeDesignEvaluated){
	lastConfigNr <- max(mergedData$CONFIG)
	lastStepNr <- mergedData$step.last
	colnames(largeDesignEvaluated)= row.names(spotConfig$alg.roi); 
	### We select only #seq.design.oldBest.size design points. These design points will be considered for re-evaluation:
	if (spotConfig$seq.design.oldBest.size <= 1) warning("spotRepeatsOcba.R: Increase spotConfig$seq.design.oldBest.size in your conf file, OCBA will not work if it is smaller than 2.");
	selection <- order(mergedData$mergedY)[1:min(nrow(mergedData$x), spotConfig$seq.design.oldBest.size)];
	ocbaData <- cbind(  mergedData$x[selection,,drop=FALSE]
					, mergedY =  as.numeric(mergedData$mergedY[selection])
					, varY = as.numeric(mergedData$varY[selection])
					, count = as.numeric(mergedData$count[selection])
					, CONFIG = as.numeric(mergedData$CONFIG[selection])
					, STEP = as.numeric(mergedData$STEP[selection])
					, SEED = as.numeric(mergedData$SEED[selection])
					) ;
	varY <- ocbaData$varY	
  
  ## check if variance >0, yields error.
  if(any(is.na(varY))|any(is.null(varY))|any(is.nan(varY))|any(varY==0)){
    stop("			
There is no variance for some point(s) in the current design.
Therefore OCBA cannot be used. Possible reasons are a target 
function without noise, or the design points are not 
evaluated repeatedly.
SPOT with OCBA makes only sense if the target function is noisy.
That means: either use spot.ocba=FALSE, or set the repeats
(init.design.repeats) to values larger than 1.

The current variance vector for the design points is: 
",paste(varY," "))
  }
  
	##
	### Based on OCBA, the budget is distributed among this subset:
	REPEATS <- spotOcba(ocbaData$mergedY, varY, ocbaData$count, spotConfig$seq.ocba.budget, iz=NA, verbose=spotConfig$io.verbosity)
	spotPrint(spotConfig$io.verbosity,1,REPEATS)
	oldD <- cbind(ocbaData, REPEATS=data.frame(REPEATS))                  
	oldD <- oldD[oldD$REPEATS>0,];
	oldD$SEED <- oldD$SEED + 1;
	oldD$mergedY <- NULL
	oldD$varY <- NULL
	oldD$count <- NULL
	### 
	additionalConfigNumbers <- nrow(largeDesignEvaluated)
	CONFIG <- lastConfigNr + 1:additionalConfigNumbers;
	REPEATS <- rep(spotConfig$init.design.repeats, additionalConfigNumbers)
	SEED <- spotConfig$alg.seed
	STEP <- lastStepNr +1
	newD <- cbind(largeDesignEvaluated
					  , CONFIG
					  , REPEATS
					  , STEP
					  , SEED)
	### Combination of old (which should be re-evaluated)  and new design 
	rbind(newD, oldD)
}

##################################################################################
#' Select repeats of design points, linearly increasing
#'
#' @param spotConfig List, containing all settings of SPOT
#' @param mergedData merged data from runs of the target function
#' @param largeDesignEvaluated sorted large design in sequential SPOT step
#' 
#' @seealso  \code{\link{spotGenerateSequentialDesign}} \code{\link{spotOcba}} \code{\link{spotRepeatsOcba}}
#' @keywords internal
###################################################################################
spotRepeats <- function(spotConfig,mergedData,largeDesignEvaluated){
	if(length(spotConfig$alg.resultColumn)>1){
		selection <- order(nds_rank(t(mergedData$mergedY)))[1:min(spotConfig$seq.design.oldBest.size,length(mergedData$CONFIG))] #TODO: this is not sufficient. hypervolume contribution, number of repeats, etc should be considered, too.
	}else{
		selection <- order(mergedData$mergedY)[1:min(spotConfig$seq.design.oldBest.size,length(mergedData$CONFIG))]
	}
	selectedData=as.data.frame(mergedData$x[selection,]) #MZ: Bugfix for 1 dimensional optimization
	colnames(selectedData)= row.names(spotConfig$alg.roi); #MZ: Bugfix for 1 dimensional optimization
	colnames(largeDesignEvaluated)= row.names(spotConfig$alg.roi); #MZ: Bugfix for 1 dimensional optimization
	oldD <- cbind(selectedData #MZ: Bugfix for 1 dimensional optimization
			, CONFIG = mergedData$CONFIG[selection]           
			, repeatsInternal = mergedData$count[selection]
			, repeatsLastConfig = mergedData$count[max(mergedData$CONFIG)] # holds the number of repeats used for the last configuration of the last step...
	)
	## new, increased number of experiments
	## definable increase function is used - see  seq.design.increase.func in spotGetOptions 
	if(!exists(spotConfig$seq.design.increase.func))stop(paste("The design increase function name", spotConfig$seq.design.increase.func, "is not found in the workspace \n
		Please make sure to load the design increase function in the workspace, or specify the correct function in spotConfig$seq.design.increase.func" ))
	totalWanted<-(eval(call(spotConfig$seq.design.increase.func, 
								max(oldD$repeatsLastConfig))));		
	## seq.design.maxRepeats is the upper bound, so the increasing repeats are limited to that maximum 
	totalWanted <- min(totalWanted,
			spotConfig$seq.design.maxRepeats, 
			na.rm=TRUE);
	## now calculate the number of repeats for those configurations, that were 
	## already evaluated before
	oldD$repeatsInternal <- totalWanted - oldD$repeatsInternal; 
	oldD = oldD[oldD$repeatsInternal>0,]	#remove those with zero repeats (i.e. reached max-repeats)

	newCONFIG <- max(mergedData$CONFIG) + 1:nrow(largeDesignEvaluated);
	newD <- cbind(  largeDesignEvaluated
			, CONFIG = newCONFIG
			, repeatsInternal = totalWanted
			, repeatsLastConfig= totalWanted);
	## if old design points have to be evaluated:
	if (sum(oldD$repeatsInternal,na.rm=TRUE) > 0){
		design <- rbind(newD,oldD)}
	## otherwise take the new design points only:
	else{
		design <- newD}
	## now replace the internal identifier with the correct one from spotConfig
	colnames(design)[colnames(design)=="repeatsInternal"] <-"REPEATS";
	## append column with current step
	design <- cbind(design,mergedData$step.last + 1);
	colnames(design)[ncol(design)] <- "STEP";
	SEED<-spotConfig$alg.seed+totalWanted-design[,"REPEATS"]
	design <- cbind(design,SEED);
	## is the following necessary?
	colnames(design)[ncol(design)] <- "SEED";
	design
}
