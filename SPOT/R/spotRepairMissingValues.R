###################################################################################
#' Repair Missing Values 
#' 
#' In case the target function that is optimized yields NA values (due to implicit constraints),
#' this function may be used to repair missing values. Currently, it can be used with Kriging
#' models only.
#'
#' @param x data frame of locations
#' @param y vector of observations
#' @param spotConfig configuration, including spotConfig$seq.na.penalty, spotConfig$seq.predictionModel.func and spotConfig$seq.modelFit
#' 
#' This function may do two things:\cr
#' 1. (optional) If seq.na.penalty is not NULL, a randomForest is trained that tries to learn the dependency of NA-occurrence on x.
#' This model will later be used to penalize potentially NA-yielding candidate solutions during surrogate-model optimization.\cr
#' 2. More importantly, this function will try to replace any NA value with the predicted mean + variance as determined
#' with the most recently used Kriging model. In case NA values occur in the initial design, they are simply removed.
#'
#' @return list with elements\cr
#' \code{x} Locations (If not repaired, locations with NA observations are removed) \cr
#' \code{y} Observations (If not repaired (due to no available model) the respective NAs are removed)\cr
#' \code{nafit} Optional fit of a randomForest model (seq.na.penalty set to positive value) that learns the occurrence of NAs for a penalty approach.
#' 
#' @seealso \code{\link{spotPenalizeMissingValues}}, \code{\link{spotRepairMissingValuesCoKriging}}
#' @export
###################################################################################
spotRepairMissingValues <- function(x, y, spotConfig){
	NAs <- is.na(y)
	fit <- NULL
	if(any(NAs)){	
		if(!is.null(spotConfig$seq.na.penalty)){ #learn the location of NA values with random forest
			fit = randomForest(x=x,y=as.factor(NAs))
		}		
		if( max(spotConfig$alg.currentResult$STEP)==0 | !is.null(spotConfig$seq.na.penalty)){ #In initial design, NA values are removed. Also, removed when they are learned with a classifier
			y <- y[which(!NAs)] 
			x <- x[which(!NAs),,drop=FALSE]
		}else{ #in later steps, NA values are replaced by predicted mean + variance
			xreplace <- x[which(NAs),,drop=FALSE]
			spotConfig$seq.model.variance <- T
			spotConfig <- eval(call(spotConfig$seq.predictionModel.func
									, NULL 
									, NULL
									, as.data.frame(xreplace)
									, spotConfig
									, spotConfig$seq.modelFit #external fit is used, model is only evaluated not build
									))
			y[which(NAs)] <- unlist(spotConfig$seq.largeDesignY + spotConfig$seq.largeDesignVar) #penalize with variance
		}
	}
	list(y=y,x=x,nafit=fit)
}

###################################################################################
#' Penalize Missing Values 
#' 
#' Values that are estimated to yield NAs will be penalized.
#'
#' @param x data frame of locations
#' @param y vector of observations
#' @param nafit randomForest fit (object), that models dependency of NAs on x
#' @param penalty (usually positive) penalty value
#'
#' @return \code{y}: Observations + penalty where applicable
#' 
#' @seealso \code{\link{spotRepairMissingValues}}
#' @export
###################################################################################
spotPenalizeMissingValues <- function(x,y,nafit,penalty){ #nafit is created by spotRepairMissingValues. penalty is spotConfig$seq.na.penalty --> which is a vector in case of MCO, else scalar
	NAs <- as.logical(predict(nafit,x))
	y[which(NAs)] <- y[which(NAs)] + penalty
	y
}

###################################################################################
#' Repair Missing Values 
#' 
#' Similar to \code{\link{spotRepairMissingValues}}, but concerns the case of
#' multi-fidelity optimization with Co-Kriging only.
#'
#' @param xe data frame of locations (expensive, fine function)
#' @param ye vector of observations (expensive, fine function)
#' @param xc data frame of locations (cheap, coarse function)
#' @param yc vector of observations (cheap, coarse function)
#' @param spotConfig configuration, including spotConfig$seq.na.penalty, spotConfig$seq.predictionModel.func and spotConfig$seq.modelFit
#' @param fitC Fit of a Kriging model (from \code{\link{forrBuilder}})
#' 
#' @return \code{y}: Observations + penalty where applicable
#' 
#' @seealso \code{\link{spotRepairMissingValues}}
#' @export
###################################################################################
spotRepairMissingValuesCoKriging <- function(xe, ye, xc, yc, spotConfig, fitC){
	le = length(ye)
	lc = length(yc)
	###
	## repair all yc value by the given fitC
	###
	NAs <- is.na(yc)
	if(any(NAs)){
		xreplace <- xc[which(NAs),,drop=FALSE]
		config <- spotConfig
		config$seq.model.variance <- T
		config$seq.modelFit <- fitC #here, only the cheap model fit can be used.
		config$seq.predictionModel.func <- "spotPredictForrester" #with the normal kriging model	
		config <- eval(call(config$seq.predictionModel.func
										, NULL 
										, NULL
										, as.data.frame(xreplace)
										, config
										, fitC #external fit is used, model is only evaluated not build
										))
		yc[which(NAs)] <- unlist(config$seq.largeDesignY + config$seq.largeDesignVar) #penalize with variance
	}
	###
	## repair NA's in ye, and move corresponding yc entries if necessary
	###
	ye1 <- ye
	listXY <- spotRepairMissingValues(xe,ye,spotConfig)
	ye = as.matrix(listXY$y)
	xe = listXY$x
	if( any(is.na(ye1)) && max(spotConfig$alg.currentResult$STEP)==0){ #In initial design, NA values are removed.
		#indf <- c(rep(F,lc-le),rep(T,length(ye)))
		xcF <- xc[1:(lc-le),,drop=F] #cheap x values which are not part of the expensive design
		xcT <- xc[(lc-le+1):lc,,drop=F]	 # "" which are part of the expensive design
		ycF <- yc[1:(lc-le),,drop=F] #same for y values
		ycT <- yc[(lc-le+1):lc,,drop=F]		
		move <- which(is.na(ye1))
		xcF <- rbind(xcF,xcT[move,,drop=F])
		ycF <- rbind(ycF,ycT[move,,drop=F])
		xcT <- xcT[-move,,drop=F]
		ycT <- ycT[-move,,drop=F]
		xc <- rbind(xcF,xcT)
		yc <- rbind(ycF,ycT)
	}	
	list(ye=ye,xe=xe,yc=yc,xc=xc)
}