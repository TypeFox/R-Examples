#
#seed    random seed (e.g. 12345)
#steps   maximum number of evolution steps (e.g. 10000)
#target  objective function threshold for priliminary evolution end (e.g. 0.0001)
#f       objective function class name (e.g. "de.fhkoeln.spot.objectivefunctions.Ball")
#n       problem dimension (e.g. 12)
#sigma0  initial step size (e.g. 1.0)
#a       step size muliplier (e.g. 1.2239)
#g       history length ( e.g. 12 == n)

###################################################################################################
#' one plus one ES function call for SPOT
#'
#' SPOT uses this function for some demos to call an one plus one evolution strategy. 
#'
#' @param spotConfig Contains the list of spot configurations, results of the algorithm can be passed to this list instead of the .res file.
#'		  spotConfig defaults to "NA", and will only be passed to the Algorithm if spotConfig$spot.fileMode=FALSE. See also: \code{\link{spotGetOptions}}
#'			Items used are: 
#'			\item{alg.currentDesign}{data frame holding the design points that will be evaluated} 
#'			\item{io.apdFileName}{name of the apd file} 
#'			\item{io.desFileName}{name of the des file} 
#'			\item{io.resFileName}{name of the res file, for logging results (if spotConfig$spot.fileMode==TRUE)}
#'			\item{spot.fileMode}{boolean, if selected with true the results will also be written to the res file, otherwise it will only be saved in the spotConfig returned by this function}
#' @return this function returns the \code{spotConfig} list with the results in spotConfig$alg.currentResult
#' @references  \code{\link{SPOT}}
###################################################################################################
spotAlgStartOnePlusOneEsJava<- function(spotConfig){
	io.apdFileName=spotConfig$io.apdFileName;
	io.desFileName=spotConfig$io.desFileName;
	io.resFileName=spotConfig$io.resFileName;	
	
	writeLines("1+1-ES run...", con=stderr());
	print(io.apdFileName)
	## read default problem design
	source(io.apdFileName,local=TRUE)
	## read doe/dace etc settings:
	if (spotConfig$spot.fileMode){ ##Check if spotConfig was passed to the algorithm, if yes the spot.fileMode is chosen with False wich means results have to be passed to spotConfig and not to res file.
		writeLines(paste("Loading design file data from::",  io.desFileName), con=stderr());
		## read doe/dace etc settings:
		des <- read.table( io.desFileName, sep=" ", header = TRUE);	
	}else{
		des <- spotConfig$alg.currentDesign; ##The if/else should not be necessary anymore, since des will always be written into the spotConfig
	}
	print(summary(des));	
	##  SIGMANULL VARA VARG REPEATS SEED
	config<-nrow(des);
	print(config);
	attach(des)	
	for (k in 1:config){
		if(des$REPEATS[k]>=1){
			for (i in 1:des$REPEATS[k]){
				##
				if (exists("SIGMANULL")){
					sigma0 <- des$SIGMANULL[k]
				}
				if (exists("VARA")){
					a <- des$VARA[k]
				}
				if (exists("VARG")){
					g <- round(des$VARG[k])
				}
				conf <- k
				if (exists("CONFIG")){
					conf <- des$CONFIG[k]
				}
				spotStep<-NA
				if (exists("STEP")){
					spotStep <- des$STEP[k]
				}
				seed <- des$SEED[k]+i			
				print(c("Config:",k ," Repeat:",i))
				callString <- paste("java -jar bin/simpleOnePlusOneES.jar", seed, steps, target, f, n, xp0, sigma0, a, g, px, py, sep = " ")
				print(callString)
				y <-system(callString, intern= TRUE)
				print(y)				
				res <- NULL
				res <- list(Y=as.numeric(as.character(y)), #converted to numeric because the java call return value is seen as a factor, which leads to problems later on.
						SIGMANULL=sigma0,
						VARA=a,
						VARG=g,
						Function=f,
						MAXITER=steps,
						DIM=n,
						TARGET=target,
						SEED=seed,
						CONFIG=conf
				)
				if (exists("STEP")){
					res=c(res,STEP=spotStep)
				} 
				res <-data.frame(res)
				if (spotConfig$spot.fileMode){ ##Log the result in the .res file, only if user didnt set fileMode==FALSE
					colNames = TRUE
					if (file.exists(io.resFileName)){
						colNames = FALSE
					}					
					## quote = false is required for JAVA
					write.table(res
						, file =  io.resFileName
						, row.names = FALSE
						, col.names = colNames
						, sep = " "              
						, append = !colNames
						, quote = FALSE
					);		
					colNames = FALSE
				}
				spotConfig$alg.currentResult=rbind(spotConfig$alg.currentResult,res);				
			} # end for i
		} # end if(des$REPEATS[k]>=1)
	}	#end for k
	detach(des)
	return(spotConfig)
}


