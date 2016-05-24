###################################################################################################
#' Interface for SANN to be tuned by SPOT
#'
#' SPOT uses this function for some demos to call the \code{\link{optim}} function with the SANN
#' method, which means Simulated Annealing. The SANN uses \code{\link{spotFuncStartBranin}} as 
#' a target function.
#' This function is needed as an interface, to ensure the right information
#' are passed from SPOT to the target algorithm(e.g. the SANN) and vice versa.
#'
#' @param spotConfig Contains the list of spot configurations, results of the algorithm can be passed to this list instead of the .res file.
#'		  spotConfig defaults to "NA", and will only be passed to the Algorithm if spotConfig$spot.fileMode=FALSE. See also: \code{\link{spotGetOptions}}
#'			Items used are: \cr \cr
#'			alg.currentDesign: data frame holding the design points that will be evaluated \cr
#'			io.apdFileName: name of the apd file \cr
#'			io.desFileName: name of the des file \cr
#'			io.resFileName: name of the res file, for logging results (if spotConfig$spot.fileMode==TRUE)\cr
#'			spot.fileMode: boolean, if selected with true the results will also be written to the res file, otherwise it will only be saved in the spotConfig returned by this function\cr
#' @return this function returns the \code{spotConfig} list with the results in spotConfig$alg.currentResult
#' @seealso  \code{\link{SPOT}} \code{\link{spot}} \code{\link{demo}} \code{\link{optim}}
#' \code{\link{spotFuncStartBranin}} \code{\link{spotAlgStartSannVar}}
#' @export
###################################################################################################
spotAlgStartSann <- function(spotConfig){
	if(exists(as.character(substitute(.Random.seed))))
		SAVESEED<-.Random.seed
	else
		SAVESEED=NULL
	io.apdFileName=spotConfig$io.apdFileName
	io.desFileName=spotConfig$io.desFileName
	io.resFileName=spotConfig$io.resFileName
	#default Values that can be changed with apd file
	func<-if(is.null(spotConfig$apd.func)){spotBraninFunction}else{spotConfig$apd.func}
  x1<-if(is.null(spotConfig$apd.x0)){c(10,10)}else{spotConfig$apd.x0}#either start point for optim, or NA for random startpoint
	maxit<-ifelse(is.null(spotConfig$apd.maxit),100,spotConfig$apd.maxit)
	parscale<-if(is.null(spotConfig$apd.parscale)){c(1,1)}else{spotConfig$apd.parscale}
	f<-ifelse(is.null(spotConfig$apd.f),"Branin",spotConfig$apd.f)
	n<-ifelse(is.null(spotConfig$apd.n),2,spotConfig$apd.n)
	lowerLimit<-if(is.null(spotConfig$apd.lower)){rep(-1,2)}else{spotConfig$apd.lower}
	upperLimit<-if(is.null(spotConfig$apd.upper)){rep(1,2)}else{spotConfig$apd.upper}
	## read problem design file
	if(file.exists(io.apdFileName)){
		source(io.apdFileName,local=TRUE)
	}
	else{
		spotWriteLines(spotConfig$io.verbosity,1,"apd File not found, defaults and conf file used")
	}	
	if (spotConfig$spot.fileMode){ 
		spotWriteLines(spotConfig$io.verbosity,1,paste("Loading design file data from::",  io.desFileName), con=stderr())
		## read doe/dace etc settings:
		des <- read.table( io.desFileName, sep=" ", header = TRUE)
	}else{
		des <- spotConfig$alg.currentDesign
	}		
	pNames <- names(des)
	config<-nrow(des)
	for (k in 1:config){
		for (i in 1:des$REPEATS[k]){
			if (is.element("TEMP", pNames)){
				temp <- des$TEMP[k]
			}
			if (is.element("TMAX", pNames)){
				tmax <- round(des$TMAX[k])
			}
			conf <- k
			if (is.element("CONFIG", pNames)){
				conf <- des$CONFIG[k]
			}
			spotStep<-NA
			if (is.element("STEP", pNames)){
				spotStep <- des$STEP[k]
			}			
			seed <- des$SEED[k]+i-1	
			set.seed(seed)
			x0<-if(is.na(x1[1])){runif(n)*(upperLimit-lowerLimit)+lowerLimit}else{x1} #random x0 if NA
			y <- optim(x0, func, method="SANN",
					control=list(maxit=maxit, temp=temp, tmax=tmax, parscale=parscale))
			res <- NULL
			res <- list(Y=y$value, TEMP=temp, TMAX=tmax, FUNCTION=f, DIM=n, SEED=seed, CONFIG=conf)
			if (is.element("STEP", pNames)){
				res=c(res,STEP=spotStep)
			} 
			res <-data.frame(res)
			if (spotConfig$spot.fileMode){ ##Log the result in the .res file, only if user didnt set fileMode==FALSE
				colNames = TRUE
				if (file.exists(io.resFileName)){
					colNames = FALSE
				}				
				write.table(res, file = io.resFileName, row.names = FALSE, 
					col.names = colNames, sep = " ", append = !colNames, quote = FALSE)
				colNames = FALSE					
			}
			spotConfig$alg.currentResult=rbind(spotConfig$alg.currentResult,res)#always log the results in spotConfig			
		}
	}	
	if(!is.null(SAVESEED))
		assign(".Random.seed", SAVESEED, envir=globalenv())
	spotConfig
}

###################################################################################################
#' Interface for SANN to be tuned robustly by SPOT
#'
#' SPOT uses this function for some demos to call the \code{\link{optim}} function with the SANN
#' method, which means Simulated Annealing. The SANN uses \code{\link{spotFuncStartBranin}} as 
#' a target function.
#' This function is needed as an interface, to ensure the right information
#' are passed from SPOT to the target algorithm(e.g. the SANN) and vice versa.
#' In contrast to \code{\link{spotAlgStartSann}} it is an interface for Pareto optimization, to optimize both the
#' performance as well as the variance of the SANN algorithm, to reach more robust results.
#'
#' @param spotConfig Contains the list of spot configurations, results of the algorithm can be passed to this list instead of the .res file.
#'		  spotConfig defaults to "NA", and will only be passed to the Algorithm if spotConfig$spot.fileMode=FALSE. See also: \code{\link{spotGetOptions}}
#'			Items used are: \cr \cr
#'			alg.currentDesign: data frame holding the design points that will be evaluated \cr
#'			io.apdFileName: name of the apd file \cr
#'			io.desFileName: name of the des file \cr
#'			io.resFileName: name of the res file, for logging results (if spotConfig$spot.fileMode==TRUE)\cr
#'			spot.fileMode: boolean, if selected with true the results will also be written to the res file, otherwise it will only be saved in the spotConfig returned by this function\cr
#' @return this function returns the \code{spotConfig} list with the results in spotConfig$alg.currentResult
#' @seealso  \code{\link{SPOT}} \code{\link{spot}} \code{\link{demo}} \code{\link{optim}}
#' \code{\link{spotFuncStartBranin}} \code{\link{spotAlgStartSann}}
#' @export
###################################################################################################
spotAlgStartSannVar <- function(spotConfig){
	if(exists(as.character(substitute(.Random.seed))))
		SAVESEED<-.Random.seed
	else
		SAVESEED=NULL
	spotConfig$alg.resultColumn=c("Y", "Ysd")
	io.apdFileName=spotConfig$io.apdFileName
	io.desFileName=spotConfig$io.desFileName
	io.resFileName=spotConfig$io.resFileName
	#default Values that can be changed with apd file
	func<-if(is.null(spotConfig$apd.func)){spotBraninFunction}else{spotConfig$apd.func}
	x1<-if(is.null(spotConfig$apd.x0)){c(10,10)}else{spotConfig$apd.x0}
	maxit<-ifelse(is.null(spotConfig$apd.maxit),100,spotConfig$apd.maxit)
	parscale<-if(is.null(spotConfig$apd.parscale)){c(1,1)}else{spotConfig$apd.parscale}
	f<-ifelse(is.null(spotConfig$apd.f),"Branin",spotConfig$apd.f)
	n<-ifelse(is.null(spotConfig$apd.n),2,spotConfig$apd.n)
	lowerLimit<-if(is.null(spotConfig$apd.lower)){rep(-1,2)}else{spotConfig$apd.lower}
	upperLimit<-if(is.null(spotConfig$apd.upper)){rep(1,2)}else{spotConfig$apd.upper}
	## read problem design file
	if(file.exists(io.apdFileName)){
		source(io.apdFileName,local=TRUE)
	}
	else{
		spotWriteLines(spotConfig$io.verbosity,1,"apd File not found, defaults used")
	}
	if (spotConfig$spot.fileMode){ 
		spotWriteLines(spotConfig$io.verbosity,1,paste("Loading design file data from::",  io.desFileName), con=stderr())
		## read doe/dace etc settings:
		des <- read.table( io.desFileName, sep=" ", header = TRUE)
	}else{
		des <- spotConfig$alg.currentDesign
	}			
	pNames <- names(des)
	config<-nrow(des)
	for (k in 1:config){
		for (i in 1:des$REPEATS[k]){
			if (is.element("TEMP", pNames)){
				temp <- des$TEMP[k]
			}
			if (is.element("TMAX", pNames)){
				tmax <- round(des$TMAX[k])
			}
			conf <- k
			if (is.element("CONFIG", pNames)){
				conf <- des$CONFIG[k]
			}
			spotStep<-NA
			if (is.element("STEP", pNames)){
				spotStep <- des$STEP[k]
			}			
			seed <- des$SEED[k]+i-1	
			set.seed(seed)
			x0<-if(is.na(x1[1])){runif(n)*(upperLimit-lowerLimit)+lowerLimit}else{x1} #random x0 if NA
			y <- optim(x0, func, method="SANN",
					control=list(maxit=maxit, temp=temp, tmax=tmax, parscale=parscale))	
			res <- list(Y=y$value, Ysd=NA, TEMP=temp, TMAX=tmax, FUNCTION=f, DIM=n, SEED=seed, CONFIG=conf)
			if (is.element("STEP", pNames)){
				res=c(res,STEP=spotStep)
			} 
			res <-data.frame(res)
			spotConfig$alg.currentResult=rbind(spotConfig$alg.currentResult,res);#always log the results in spotConfig			
		}
		ysd<-sd(spotConfig$alg.currentResult[which(spotConfig$alg.currentResult$CONFIG==des$CONFIG[k]),]$Y)
		spotConfig$alg.currentResult[which(spotConfig$alg.currentResult$CONFIG==des$CONFIG[k]),]$Ysd=ysd
	}
	#In this function, logging can not be done in between, since variance is only known in the end
	if (spotConfig$spot.fileMode){ ##Log the result in the .res file, only if user didnt set fileMode==FALSE
		colNames = TRUE
		if (file.exists(io.resFileName)){
			colNames = FALSE
		}				
		write.table(spotConfig$alg.currentResult[which(spotConfig$alg.currentResult$STEP==spotStep),],
					file = io.resFileName, row.names = FALSE, 
					col.names = colNames, sep = " ", append = !colNames, quote = FALSE)
		colNames = FALSE					
	}
	if(!is.null(SAVESEED))
		assign(".Random.seed", SAVESEED, envir=globalenv())
	spotConfig
}