
###################################################################################################
#' Branin function call interface for SPOT, deprecated
#'
#' SPOT uses this function for some demos to call the \code{\link{spotBraninFunction}} function 
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
#' \code{\link{spotBraninFunction}}
#' @keywords internal
###################################################################################################
spotFuncStartBranin <- function(spotConfig){
	pdFile=spotConfig$io.apdFileName
	resFileName=spotConfig$io.resFileName
	desFileName=spotConfig$io.desFileName
	if(is.null(spotConfig$spot.noise)){spotConfig$spot.noise=10.0}
	if(is.null(spotConfig$spot.noise.type)){spotConfig$spot.noise.type="weighted"}
	if(is.null(spotConfig$spot.noise.minimum.at.value)){spotConfig$spot.noise.minimum.at.value=0.0}
	
	if (spotConfig$spot.fileMode){ 
		spotWriteLines(spotConfig$io.verbosity,1,paste("Loading design file data from::",  desFileName), con=stderr())
		## read doe file
		des <- read.table(desFileName, sep=" ", header = TRUE)
	}else{
		des <- spotConfig$alg.currentDesign 
	}
	spotPrint(spotConfig$io.verbosity,1,summary(des))
	spotWriteLines(spotConfig$io.verbosity,1,"Branin run...", con=stderr())
	spotPrint(spotConfig$io.verbosity,1,pdFile)
	#default Values that can be changed with apd file
	noise<-spotConfig$spot.noise
	noise.type <- spotConfig$spot.noise.type
	spot.noise.minimum.at.value <- spotConfig$spot.noise.minimum.at.value
	f<-"Branin"
	n<-2
	## read problem design file
	if(file.exists(pdFile)){
		source(pdFile,local=TRUE)
	}
	##  VARX1 VARX2 REPEATS SEED
	config<-nrow(des)
	spotPrint(spotConfig$io.verbosity,1,config)
	for (k in 1:config){
		for (i in 1:des$REPEATS[k]){
			##
			if (!is.null(des$VARX1)){
				x1 <- des$VARX1[k]
			}
			if (!is.null(des$VARX2)){
				x2 <- des$VARX2[k]
			}
			conf <- k
			if (!is.null(des$CONFIG)){
				conf <- des$CONFIG[k]
			}
			if (!is.null(des$STEP)){
				step <- des$STEP[k]
			}
			seed <- des$SEED[k]+i-1			
			spotPrint(spotConfig$io.verbosity,1,c("Config:",k ," Repeat:",i))
			y <- spotBraninFunction(c(x1,x2))
			## add noise
			y <- y + spotCalcNoise(y, noise=noise, noise.type=noise.type, spot.noise.minimum.at.value=spot.noise.minimum.at.value);
			spotPrint(spotConfig$io.verbosity,1,y)
			res <- NULL
			res <- list(Y=y,					
					VARX1=x1,
					VARX2=x2,
					Function=f,					
					DIM=n,
					STEP=step,
					SEED=seed,
					CONFIG=conf
			)
			res <-data.frame(res)
			if (spotConfig$spot.fileMode){ ##Log the result in the .res file, only if user didnt set fileMode==FALSE
				colNames = TRUE
				if (file.exists(resFileName)){
					colNames = FALSE
				}				
				## quote = false is required for JAVA
				write.table(res
						, file = resFileName
						, row.names = FALSE
						, col.names = colNames
						, sep = " "              
						, append = !colNames
						, quote = FALSE
				);		
			}
			spotConfig$alg.currentResult=rbind(spotConfig$alg.currentResult,res)	#always log the results in spotConfig				
		}			
	}	
	spotConfig
}

###################################################################################################
#' Stairlike Branin function call interface for SPOT, deprecated
#'
#' SPOT uses this function for some demos to call the \code{\link{spotBraninFunction}} function 
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
#' \code{\link{spotBraninFunction}}
#' @keywords internal
###################################################################################################
spotFuncStartDisBranin <- function(spotConfig){
	pdFile=spotConfig$io.apdFileName
	resFileName=spotConfig$io.resFileName
	desFileName=spotConfig$io.desFileName
	if (spotConfig$spot.fileMode){
		spotWriteLines(spotConfig$io.verbosity,1,paste("Loading design file data from::",  desFileName), con=stderr())
		## read doe/dace etc settings:
		des <- read.table(desFileName, sep=" ", header = TRUE)
	}else{
		des <- spotConfig$alg.currentDesign
	}
	spotPrint(spotConfig$io.verbosity,1,summary(des))
	spotWriteLines(spotConfig$io.verbosity,1,"Branin run...", con=stderr())
	spotPrint(spotConfig$io.verbosity,1,pdFile)
	#default Values that can be changed with apd file
	noise<-spotConfig$spot.noise
	noise.type <- spotConfig$spot.noise.type
	spot.noise.minimum.at.value <- spotConfig$spot.noise.minimum.at.value
	f<-"Branin"
	n<-2;
	## read problem design file
	if(file.exists(pdFile)){
		source(pdFile,local=TRUE)
	}
	##  VARX1 VARX2 REPEATS SEED
	config<-nrow(des);
	spotPrint(spotConfig$io.verbosity,1,config);
	for (k in 1:config){
		for (i in 1:des$REPEATS[k]){
			##
			if (!is.null(des$VARX1)){
				x1 <- des$VARX1[k]
			}
			if (!is.null(des$VARX2)){
				x2 <- des$VARX2[k]
			}
			conf <- k
			if (!is.null(des$CONFIG)){
				conf <- des$CONFIG[k]
			}
			if (!is.null(des$STEP)){
				step <- des$STEP[k]
			}
			seed <- des$SEED[k]+i-1			
			spotPrint(spotConfig$io.verbosity,1,c("Config:",k ," Repeat:",i))
			y <- spotBraninFunction(c(x1,x2))
			if(y>150){y=150}
			else if(y<=150&y>140){y=160.0}
			else if(y<=140&y>130){y=110.1}
			else if(y<=130&y>110){y=50.2}
			else if(y<=130&y>110){y=33.3}
			else if(y<=110&y>101){y=38.4}
			else if(y<=101&y>88){y=66.7}
			else if(y<=88&y>50){y=47.3}
			else if(y<=50&y>45){y=33.1}
			else if(y<=45&y>32){y=25.9}
			else if(y<=32&y>25){y=28.6}
			else if(y<=25&y>18){y=21.3}
			else if(y<=18&y>12){y=24.5}
			else if(y<=12&y>8){y=21.2}
			else if(y<=8&y>5){y=10.1}
			else if(y<=5&y>0){y=5}
			else if(y<=0&y>0.2){y=4}
			else if(y<=0.2){y=0.4}
			## add noise
			y <- y + spotCalcNoise(y, noise=noise, noise.type=noise.type, spot.noise.minimum.at.value=spot.noise.minimum.at.value);
			spotPrint(spotConfig$io.verbosity,1,y)
			res <- NULL
			res <- list(Y=y,					
					VARX1=x1,
					VARX2=x2,
					Function=f,					
					DIM=n,
					STEP=step,
					SEED=seed,
					CONFIG=conf
			)
			res <-data.frame(res)
			if (spotConfig$spot.fileMode){ ##Log the result in the .res file, only if user didnt set fileMode==FALSE
				colNames = TRUE
				if (file.exists(resFileName)){
					colNames = FALSE
				}				
				## quote = false is required for JAVA
				write.table(res
						, file = resFileName
						, row.names = FALSE
						, col.names = colNames
						, sep = " "              
						, append = !colNames
						, quote = FALSE)		
			}
			spotConfig$alg.currentResult=rbind(spotConfig$alg.currentResult,res);	#always log the results in spotConfig				
		}			
	}	
	spotConfig
}
