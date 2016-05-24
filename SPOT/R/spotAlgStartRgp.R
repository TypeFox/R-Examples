###################################################################################################
#' Target function for RGP / SPOT tuning
#'
#' This target function is used by \code{\link{spotAlgStartRgp}} to start rgp with different population or 
#' tournament sizes.
#'
#' @param populationSize the \code{populationSize} parameter, as used by \code{symbolicRegression()} in \code{rgp}
#' @param tournamentSize the \code{tournamentSize} parameter, as used by \code{makeTournamentSelection()} in \code{rgp}
#' @param time the time for the stopping criterion of the rgp run, in seconds
#'
#' @return this function returns the RMSE (fitness) of the best individual in the population
#'
#' @seealso \code{\link{spotAlgStartRgp}}
#' @export
#' @keywords internal
###################################################################################################
spotRgpTargetFunction <- function(populationSize = 100, tournamentSize = 10, time=120) {
  data1 <- {
    x1 <- seq(0, 4*pi, length.out=201)
    x2 <- seq(0, 4*pi, length.out=201)
    y <- sin(x1) + cos(2*x2)
    data.frame(y=y, x1=x1, x2=x2)
  }

  mdl <- rgp::symbolicRegression(y ~ x1 + x2,
                            data = data1,
                            populationSize = populationSize,
                            #selectionFunction = makeTournamentSelection(tournamentSize = tournamentSize), #deprecated
							searchHeuristic = rgp::makeArchiveBasedParetoTournamentSearchHeuristic(popTournamentSize = tournamentSize),
                            functionSet = rgp::arithmeticFunctionSet,
                            stopCondition = rgp::makeTimeStopCondition(time))

  ## Calculate RMSE (fitness) of best individual in population
  bestFitness <- min(sapply(mdl$population, mdl$fitnessFunction))

  return (bestFitness)
}
###################################################################################################
#' Interface for RGP to be tuned by SPOT
#'
#' SPOT uses this function for some demos to call the \code{symbolicRegression} function
#' from the \code{rgp} package.
#' This function is needed as an interface, to ensure the right information
#' are passed from SPOT to the target algorithm (i.e. RGP) and vice versa.
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
#' @seealso  \code{\link{SPOT}} \code{\link{spot}} \code{\link{demo}}
#' @export
###################################################################################################
spotAlgStartRgp <- function(spotConfig){
	if(exists(as.character(substitute(.Random.seed))))
		SAVESEED<-.Random.seed
	else
		SAVESEED=NULL
	spotInstAndLoadPackages("rgp") #installs and requires() rgp
	io.apdFileName=spotConfig$io.apdFileName
	io.desFileName=spotConfig$io.desFileName
	io.resFileName=spotConfig$io.resFileName	
	#default Values that can be changed with apd file
	populationSize <- 100
	tournamentSize <- 10
	time<-120
	f <- "sincos2d"; ## TODO set name of test function here!
	## read problem design file
	if(file.exists(io.apdFileName)) {
		source(io.apdFileName,local=TRUE)
	}
	else{
		spotWriteLines(spotConfig$io.verbosity,1,"apd File not found, defaults used");
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
			if (is.element("POPULATIONSIZE", pNames)){
				populationSize <- des$POPULATIONSIZE[k]
			}
			if (is.element("TOURNAMENTSIZE", pNames)){
				tournamentSize <- round(des$TOURNAMENTSIZE[k])
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
      ## force valid parameterization (repair)...
      quantize <- function(n, q) n - (n %% q)
      populationSize <- quantize(as.integer(populationSize), 4L) ## populationSize should be a multiple of 4
      populationSize <- if (populationSize <= 1000L) populationSize else 1000L
      populationSize <- if (populationSize >= 20L) populationSize else 20L
      tournamentSize <- quantize(as.integer(tournamentSize), 4L) ## tournamentSize must be a multiple of 4
      tournamentSize <- if (tournamentSize <= populationSize) tournamentSize else populationSize
      tournamentSize <- if (tournamentSize >= 4L) tournamentSize else 4L
			y <- spotRgpTargetFunction(populationSize, tournamentSize,time)
			res <- NULL
			res <- list(Y=y, POPULATIONSIZE=populationSize, TOURNAMENTSIZE=tournamentSize, FUNCTION=f, SEED=seed, CONFIG=conf)
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
