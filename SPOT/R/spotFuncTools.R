###################################################################################
#' Calculate Noise
#'
#' Calculates Noise for a given Function-Value (y) dependant on the
#' noise intensity (noise) and the way of calculating noise (noise.type)
#'
#' @param y the function value where noise will be added to
#' @param noise noise magnitude (absolute or weighted)
#' @param noise.type type of noise can either be "weighted" or "constant"
#' @param spot.noise.minimum.at.value Voice magnitude at minimum
#' 
#' @return numeric \cr
#' holding the noise Value
#' @export
#' @keywords internal
####################################################################################
spotCalcNoise <- function(y, noise=0.0, noise.type="weighted", spot.noise.minimum.at.value=0.0){
	
	noiseValue <- 0	
	if (noise == 0)return (noiseValue)

	if (noise.type=="weighted"){
		noiseValue <- (y - spot.noise.minimum.at.value)*noise*rnorm(1)/100
	}else if (noise.type=="constant"){
		noiseValue <- noise*rnorm(1)
	}
	
	noiseValue
}

###################################################################################
#' Write Result File
#'
#' Writes a new row to the result file and adds it to the existing data.frame
#'
#' @param res result to write
#' @param file name of result file
#' @param mode only write to file if mode is TRUE
#' @param cres res will be appended to this existing data frame
#' 
#' @return numeric \cr
#' holding the noise Value
#' @export
#' @keywords internal
####################################################################################
spotWriteRes <- function(res, file, mode, cres){
	if (mode){ ##Log the result in the .res file, only if user did not set fileMode==FALSE
		colNames = TRUE
		if (file.exists(file)){
			colNames = FALSE
		}				
		## quote = false is required for JAVA
		write.table(res
				, file = file
				, row.names = FALSE
				, col.names = colNames
				, sep = " "    
				, append = !colNames
				, quote = FALSE)	
	}
	rbind(cres,res)	
}

###################################################################################
#' Write Result File
#'
#' Writes a new row to the result file and adds it to the existing data.frame
#'
#' @param y result value passed from target function
#' @param colNames column names of the results in the target data frame
#' @param res target data frame
#' 
#' @return list \cr
#' containing the extended res data frame and the corrected colNames
#' @export
#' @keywords internal
####################################################################################
spotPrepareResult <- function (y,colNames,res){
	if(length(y)>1 && length(colNames)==1){
		colNames <-  paste(colNames,".", 1:length(y), sep="")
	}
	res[colNames] <- unname(y)
	list(res=res,colNames=colNames)
}	

###################################################################################################
#' Interface for Target Functions
#'
#' SPOT uses this function to call functions passed to \code{\link{spotOptim}} or \code{\link{spot}} like they would be passed to optim().
#' That means, it will be used whenever an actual function is passed instead of a string. When a string is passed
#' the string itself will contain the interface to use.
#' This function is needed as an interface, to ensure the right information
#' are passed from SPOT to the target function. It can handle single and multi criteria target functions, e.g. functions that return numerics or vectors of numerics.
#'
#' @param spotConfig Contains the list of spot configurations, results of the algorithm can be passed to this list instead of the .res file.
#'		  spotConfig defaults to "NA", and will only be passed to the Algorithm if spotConfig$spot.fileMode=FALSE. See also: \code{\link{spotGetOptions}}
#'			Items used are: \cr \cr
#'			alg.currentDesign: data frame holding the design points that will be evaluated \cr
#'			io.apdFileName: name of the apd file \cr
#'			io.desFileName: name of the des file \cr
#'			io.resFileName: name of the res file, for logging results (if spotConfig$spot.fileMode==TRUE)\cr
#'			spot.fileMode: boolean, if selected with true the results will also be written to the res file, otherwise it will only be saved in the spotConfig returned by this function\cr
#'			spotConfig$alg.tar.func target function of type y=f(x,...)
#' @param ... additional parameters to be passed on to target function: spotConfig$alg.tar.func
#' @return this function returns the \code{spotConfig} list with the results in spotConfig$alg.currentResult
#' @seealso  \code{\link{SPOT}} \code{\link{spot}} \code{\link{demo}} \code{\link{optim}}
#' \code{\link{spotOptim}}
#' @export
###################################################################################################
#TODO: Frage: Wie sollen initial defaults gehandelt werden. Voorallem: Noise, OCBA, init/seq Repeats.
spotOptimInterface <- function(spotConfig,...){
	if(exists(as.character(substitute(.Random.seed))))
		SAVESEED<-.Random.seed
	else
		SAVESEED=NULL
	
	if (spotConfig$spot.fileMode){ 
		spotWriteLines(spotConfig$io.verbosity,1,paste("Loading design file data from::",  spotConfig$io.desFileName), con=stderr())
		## read doe/dace etc settings:
		des <- read.table(spotConfig$io.desFileName, sep=" ", header = TRUE)
	}else{
		des <- spotConfig$alg.currentDesign; 
	}
	if(is.null(spotConfig$seq.co.included))spotConfig$seq.co.included=FALSE
	spotPrint(spotConfig$io.verbosity,1,summary(des))
	spotWriteLines(spotConfig$io.verbosity,1,"spotOptimInterface...", con=stderr())
	## read problem design file   #TODO: this has no effect. values from apd are not used.
	#if(file.exists(spotConfig$io.apdFileName)){
	#	source(spotConfig$io.apdFileName,local=TRUE)
	#}
	config<-nrow(des)
	pNames <- rownames(spotConfig$alg.roi)
	ndim <- nrow(spotConfig$alg.roi)
	spotPrint(spotConfig$io.verbosity,1,config)
	for (k in 1:config){
		for (i in 1:des$REPEATS[k]){			
			x <- as.matrix(des[k,,drop=FALSE][pNames])
			conf <- k
			if (!is.null(des$CONFIG)){
				conf <- des$CONFIG[k]
			}
			if(!is.na(spotConfig$alg.seed)){ #only use seed if seed is desired (not for deterministic target functions)
				seed <- des$SEED[k]+i-1	
				set.seed(seed)	
			}
			else{seed=NA}
			spotPrint(spotConfig$io.verbosity,1,c("Config:",k ," Repeat:",i))
			if(!is.function(spotConfig$alg.tar.func)){stop("spotConfig$alg.tar.func is not a function. \n Please specify a function for this variable if you use spotOptimInterface.\n Else use your own custom interface")}
			res <-  spotConfig$alg.tar.func(x,...)#spotConfig$alg.func.tar(x)#, noise=noise, noise.type=noise.type, spot.noise.minimum.at.value=spot.noise.minimum.at.value)
			if(spotConfig$seq.co.included){ #  if expensive evaluations include the evaluation of the cheap design
				y=res$y
				yco=res$yco
			}else{
				y=res
			}			
			spotPrint(spotConfig$io.verbosity,1,y)
			res <- data.frame(XDIM=ndim,YDIM=length(y),STEP=des$STEP[k],SEED=seed,CONFIG=conf,TIME=gsub("\\s","_",Sys.time()))
			res[pNames]=x
			tmp <- spotPrepareResult(y,spotConfig$alg.resultColumn,res)
			res <- tmp$res
			spotConfig$alg.resultColumn <- tmp$colNames
			spotConfig$alg.currentResult <- spotWriteRes(res,spotConfig$io.resFileName,spotConfig$spot.fileMode,spotConfig$alg.currentResult)
					
			#co-results at evaluated positions
			if(is.function(spotConfig$seq.co.func)){
				if(spotConfig$seq.co.included){ #  if expensive evaluations include the evaluation of the cheap design
					y <- yco
				}else{				
					y <-  spotConfig$seq.co.func(x,...)#spotConfig$alg.func.tar(x)#, noise=noise, noise.type=noise.type, spot.noise.minimum.at.value=spot.noise.minimum.at.value)
				}
				res$indf=T #indf: indicate that this point in the co-result is also available as fine function evaluation
				tmp <- spotPrepareResult(y,spotConfig$alg.resultColumn,res)
				res <- tmp$res
				spotConfig$alg.resultColumn <- tmp$colNames			
				spotConfig$alg.currentCoResult <- spotWriteRes(res,spotConfig$io.coResFileName,spotConfig$spot.fileMode,spotConfig$alg.currentCoResult)					
			}					
		}			
	}	
	if(is.function(spotConfig$seq.co.func)){ #evaluation of a co-design-function. #TODO: have a set of nested co design functions
		if(is.null(spotConfig$seq.co.design.size))spotConfig$seq.co.design.size=100
		if(nrow(spotConfig$alg.currentResult)==nrow(spotConfig$alg.currentCoResult) & spotConfig$seq.co.design.size > 0){ #only if not available already			
			if(is.null(spotConfig$seq.co.design.func))spotConfig$seq.co.design.func="spotCreateDesignLhd"
			if(is.null(spotConfig$seq.co.design.retries))spotConfig$seq.co.design.retries=100
			if(is.null(spotConfig$seq.co.design.repeats))spotConfig$seq.co.design.repeats=1
			if(is.null(spotConfig$seq.co.infill))spotConfig$seq.co.infill=5; #number of points created by infill criterion for co function #todo
			if(is.null(spotConfig$seq.co.recalculate))spotConfig$seq.forr.co.recalculate=FALSE
			
			#cheap design:
			tmpConfig <- spotConfig
			#tempdes<-des
			tmpConfig$nested.design <- as.matrix(des[pNames])
			des  <- (eval(call(spotConfig$seq.co.design.func, 
												tmpConfig, 
												spotConfig$seq.co.design.size, 
												spotConfig$seq.co.design.retries)));
												
			#design should have columns: xNames, REPEATS, CONFIG, STEP, (SEED)
			des$STEP=rep(0,spotConfig$seq.co.design.size)
			des$CONFIG=1:spotConfig$seq.co.design.size
			des$REPEATS=rep(spotConfig$seq.co.design.repeats,spotConfig$seq.co.design.size)
			config<-nrow(des)
			#spotConfig$seq.cres <- spotCallCoFunction(spotConfig,NULL,xc1,spotConfig$seq.co.func)		
			#
			################################
			#now the additional evaluations on low fidelty (BMP) function only
			#
			for (k in 1:config){
				for (i in 1:des$REPEATS[k]){
					x <- as.matrix(des[k,,drop=FALSE][pNames])
					conf <- k
					if (!is.null(des$CONFIG)){
						conf <- des$CONFIG[k]
					}
					if(!is.na(spotConfig$alg.seed)){ #only use seed if seed is desired (not for deterministic target functions)
						seed <- spotConfig$alg.seed+i-1	
						set.seed(seed)	
					}
					else{seed=NA}
					y <-  spotConfig$seq.co.func(x,...)
					res <- data.frame(XDIM=ndim,YDIM=length(y),STEP=des$STEP[k],SEED=seed,CONFIG=conf,TIME=gsub("\\s","_",Sys.time()),indf=F) #indf: indicate that this point in the co-result is not available as fine function evaluation
					res[pNames]=x
					tmp <- spotPrepareResult(y,spotConfig$alg.resultColumn,res)
					res <- tmp$res
					spotConfig$alg.resultColumn <- tmp$colNames					
					spotConfig$alg.currentCoResult <- spotWriteRes(res,spotConfig$io.coResFileName,spotConfig$spot.fileMode,spotConfig$alg.currentCoResult)	
				}		
			}			
		}	
	}	
	if(!is.null(SAVESEED))
		assign(".Random.seed", SAVESEED, envir=globalenv())
	spotConfig
}
