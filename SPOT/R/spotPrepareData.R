###################################################################################
#' Get Raw Result Data 
#'
#' spotGetRawResData is based on \code{\link{spotPrepareData}}  
#'
#' The result (.res) file is read and the data is returned as is 
#' 
#' @param spotConfig the list of all parameters is given, used ones are: 
#' 		\code{spotConfig$io.resFileName:} result file where data are read from
#' 		\code{spotConfig$io.columnSep:} column separator
#' parameter to the function call \code{\link{spotPrepareData}} 
#' 
#'
#' @return data.frame \code{rawResData}  \cr
#' - \code{rawResData} contains values from the result file 
#'
#' @seealso \code{\link{spotPrepareData}} 
#' @export
#' @keywords internal
####################################################################################
spotGetRawResData<- function(spotConfig){
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotGetRawResData")
	## Load .res-file data, the result of the alg run
	spotWriteLines(spotConfig$io.verbosity
					, 2
					, paste("Loading result data from::", spotConfig$io.resFileName)
					, con=stderr())	
	rawResData <- read.table(spotConfig$io.resFileName
			, sep=spotConfig$io.columnSep
			, header = TRUE	
			, stringsAsFactors = TRUE)
	if (length(spotConfig$alg.resultColumn)==1){ #If result column setting is one dimensional, check if more columns starting with that name exist, like Y.1, Y.2, etc.
		tmp<-length(grep(paste(spotConfig$alg.resultColumn,".",sep=""),names(rawResData)))
		if(tmp>1){
			spotConfig$alg.resultColumn<- names(rawResData)[grep(paste(spotConfig$alg.resultColumn,".",sep=""),names(rawResData))]
		}
	}		
	spotWriteLines(spotConfig$io.verbosity,2,"  Leaving spotGetRawResData")
	list(conf=spotConfig,rawD=rawResData)
}

###################################################################################
#' Get Raw Data Matrix B  
#'
#' The result (.res) file is read and the data are prepared for further processing
#' results in "Matrix B" that consists of the y-vector of the results and the x-Matrix
#' of inputs that are bound together to a matrix B
#' 
#' @param spotConfig the list of all parameters is given, also used as input
#' parameter to the function call spotGetRawResData(spotConfig) the other ones in use are: 
#' 		\code{spotConfig$alg.roi:} the roi data provided as matrix
#' 		\code{spotConfig$alg.resultColumn:} name of the result column 
#'
#' @return data.frame \code{B}  \cr
#' - \code{B} holds a column "y" with the results
#' and all columns with column-names derived from .roi file (should be the parameters
#' of the algorithm. Values are sorted with respect to the y values (increasing) 
#'
#' @seealso \code{\link{spotPrepareData}} 
#' @export
#' @keywords internal
####################################################################################
spotGetRawDataMatrixB <- function(spotConfig){
	## read data frame from res file
	if(spotConfig$spot.fileMode){
		res<- spotGetRawResData(spotConfig)
		spotConfig<-res$conf
		rawData<-res$rawD
		res<-NULL
	}else{
		rawData=spotConfig$alg.currentResult 
	}   
	y <- rawData[,spotConfig$alg.resultColumn,drop=FALSE]
	A <- rawData
	if(length(spotConfig$alg.resultColumn) > 1){B <- data.frame(A[order(y[,1],y[,2],decreasing=FALSE),])}
	else{B <- data.frame(A[order(y,decreasing=FALSE),])}
	B
}

###################################################################################
#' Get Merged Data Matrix B
#'
#' The merged data that are the result of \code{\link{spotPrepareData}}must be the first input parameter
#' These data prepared and results in "Matrix B" that consists of the y-vector 
#' of the results and the x-Matrix of inputs that are bound together to a matrix B
#' 
#' @param mergedData the Data prepared as done with \code{\link{spotPrepareData}}
#' @param spotConfig the list of all parameters is given, but used only as input
#' parameter to the function call of code{\link{spotPrepareData}} 
#' 
#' @return data.frame \code{B}  \cr
#' - \code{B} holds a column "y" with the results
#' and all columns with column-names derived from .roi file (should be the parameters
#' of the algorithm. Values are sorted with respect to the y values (increasing) 
#'
#' @seealso \code{\link{spotPrepareData}} 
#' @export
#' @keywords internal
####################################################################################
spotGetMergedDataMatrixB <- function(mergedData, spotConfig){
	## Note: requires pre-processing via spotPrepareData()
	##       does not work on raw data	
	## extract parameter names
	pNames <- row.names(spotConfig$alg.roi)
	y <- mergedData$mergedY
	## data frame of parameter values	
	x <- as.matrix(mergedData$x)
	A <- cbind(y,x)        
	colnames(A)<- c(spotConfig$alg.resultColumn,pNames)
	if(!is.null(dim(y))){B <- data.frame(A[order(y[,1],y[,2],decreasing=FALSE),,drop=FALSE])}
	else{B <- data.frame(A[order(y,decreasing=FALSE),,drop=FALSE])}
	B
}

###################################################################################
#' Prepare Data for SPOT
#' 
#' spotPrepareData prepares the data from .res-file  
#'
#' The result (.res) file is read and the data are prepared for further processing
#' results in an unsorted list of several values
#' 
#' @param spotConfig the list of all parameters is given, used parameters are: 
#' 		\code{spotConfig$io.resFileName}
#' 		\code{spotConfig$io.columnSep}
#' 		\code{spotConfig$alg.roi}
#' 		\code{spotConfig$alg.resultColumn}
#' 		\code{spotConfig$seq.transformation.func()}
#' 		\code{spotConfig$seq.merge.func()}
#'
#' @return list \code{resultList}  \cr
#' - \code{resultList} is a list of values (not sorted in increasing order!)
#' 		\code{x}: vector of vector (=matrix) of input values, the "points"
#' 		\code{mergedY}: vector of merged fitness values, the result column (each row of x has a corresponding Y-value)
#'		\code{count}: the number of repeats
#' 		\code{CONFIG}: the unique CONFIG number for each unique design point (sdev=NA if count=1)
#' 		\code{pNames}: names of the parameters as derived from roi-file
#' 		\code{step.last}: the maximum number of the step-column of the result-file
#' 		\code{STEP}: SPOT STEP when design point was created  
#'
#' @seealso \code{\link{SPOT}} \code{\link{spot}}
#' @export
#' @keywords internal
####################################################################################
spotPrepareData <- function(spotConfig){
	spotWriteLines(spotConfig$io.verbosity,2,"  Entering spotPrepareData")
	if(spotConfig$spot.fileMode){
		res<- spotGetRawResData(spotConfig)
		spotConfig<-res$conf
		rawData<-res$rawD
		spotConfig$alg.currentResult<-rawData
		res<-NULL
	}else{
	  rawData=spotConfig$alg.currentResult;
	}                                   # WK: if-else needed if a call 'spot(confFileName,"rep")' shall succeed (!)
	if (!any(names(rawData)=="CONFIG"))
  	   stop("Error: Result file is missing the required column CONFIG!")
    ## extract parameter names
    pNames <- row.names(spotConfig$alg.roi);
	step.last <- NA;
	if (any(names(rawData)=="STEP"))
   		step.last <- max(rawData["STEP"])
	else{
		warning("Warning: user algorithm does not write STEP column to result file, 0-vector added")
		step.last=0
		rawData$STEP<-rep(0, length(rawData$CONFIG)) 
	}
	z <- split(rawData[,spotConfig$alg.resultColumn], rawData$CONFIG)
	
	if(length(spotConfig$alg.resultColumn)==1)
	{
		fs <- function(xs) { spotConfig$seq.transformation.func(spotConfig$seq.merge.func(xs)) }
		mergedY <- sapply(z,fs)
        varY <- sapply(z,var)
		count <- sapply(z,length)
	}
	else{
		fs <- function(xs) { spotConfig$seq.transformation.func(apply(xs,2,spotConfig$seq.merge.func))}
		fnvar<-function(x){apply(x,2,var)}
		varY <- sapply(z,fnvar)#sapply(as.data.frame(z),var);
		mergedY <- t(sapply(z,fs))
		count <- sapply(z,dim)[1,]
	}
	mergedCONFIG <- sapply(split(rawData$CONFIG, rawData$CONFIG),min)
	mergedSTEP <- sapply(split(rawData$STEP, rawData$CONFIG),min)
	### added Seed 6.Jan 2011: x
        mergedSEED <- sapply(split(rawData$SEED, rawData$CONFIG),max) #MZBUGFIX: max instead of min. maximum used seed is incremented. minimum creates repeated evaluation of same seed (at least in ocba)
	#x1 <- as.data.frame(t(sapply(split(rawData[,pNames], rawData$CONFIG),mean)));	# This expression results into errors if just one variable is in the ROI
	if (length(pNames)==1){ 
		#x <- as.data.frame(sapply(split(rawData[,pNames], rawData$CONFIG),mean)) #mean on dataframes gives warnings in R-14.x
		x<- as.data.frame(unique(cbind(rawData[,pNames], rawData$CONFIG))[,1])
		names(x)<-pNames
	}
	else{
		#browser()
		x<-unique(cbind(rawData[,pNames], rawData$CONFIG))[,pNames]
		#x <- as.data.frame(t(sapply(split(rawData[,pNames], rawData$CONFIG),mean))) #mean on dataframes gives warnings in R-14.x
	}
	resultList<-list( x = x
         , mergedY = mergedY
         , varY = varY                
         , count = count          
         , CONFIG = mergedCONFIG 
         , pNames = pNames
         , step.last = step.last
         , STEP = mergedSTEP
         , SEED = mergedSEED
         )
	spotWriteLines(spotConfig$io.verbosity,2,"  Leaving spotPrepareData")
	resultList
}


###################################################################################
#' Prepare Data As Matrix C
#'
#' The result (.res) file is read and the data are prepared for further processing
#' results in "Matrix C" that consists of the y-vector of the results and the x-Matrix
#' of inputs plus the columns "count", "sdev" and "CONFIG", that are bound together 
#' to a matrix C
#' 
#' @param spotConfig the list of all parameters is given, but used only as input
#' parameter to the function call of \code{\link{spotPrepareData}}  
#'
#' @return Matrix \code{C} \cr
#' - \code{C} holds a column "y" with the results
#' and all columns with column-names derived from .roi file (should be the parameters
#' of the algorithm) plus the columns "count", "sdev" and "CONFIG"
#
#' @seealso \code{\link{spotPrepareData}} 
#' @keywords internal
####################################################################################
spotPrepareDataAsMatrixC <- function(spotConfig){ #TODO this function is used NOWHERE ?
	algResults<-spotPrepareData(spotConfig)
	x <- as.matrix(algResults$x)
	y <- algResults$Y
	A <- cbind(y,x,count=algResults$count,CONFIG=algResults$CONFIG)       
	data.frame(A[order(y,decreasing=FALSE),])	
}
