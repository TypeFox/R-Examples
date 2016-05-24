## Experimental research in evolutionary computation
## author: thomas.bartz-beielstein@fh-koeln.de
## http://www.springer.com/3-540-32026-1
##
## Copyright (C) 2003-2010 T. Bartz-Beielstein and C. Lasarczyk
## This program is free software;
## you can redistribute it and/or modify it under the terms of the 
## GNU General Public License as published by the Free Software Foundation; 
## either version 3 of the License,
## or (at your option) any later version.
## This program is distributed in the hope that it will be useful, 
## but WITHOUT ANY WARRANTY; without even the implied warranty of 
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## See the GNU General Public License for more details.
## You should have received a copy of the GNU General Public License along
##  with this program; if not, see <http://www.gnu.org/licenses/>.
##

# Package Description for Roxygen:
#' Sequential Parameter Optimization Toolbox in R
#'
#' SPOT is a package for R, using statistic models to find
#' optimal parameters for optimization algorithms. SPOT is a very flexible and 
#' user oriented tool box for parameter optimization. The flexibility has its 
#' price: to fully use all the possibilities of flexibility the user is requested 
#' to look at a large number of spot-parameters to change. The good news is, that 
#' some defaults are given that already work perfectly well for 90 percent of the users.
#'
#' \tabular{ll}{
#' Package: \tab SPOT\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.5543\cr
#' Date: \tab 2015-04-24\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @name SPOT-package
#' @aliases SPOT
#' @docType package
#' @title Sequential Parameter Optimization Toolbox in R
#' @author Thomas Bartz-Beielstein \email{thomas.bartz-beielstein@@fh-koeln.de} with contributions from: J. Ziegenhirt, W.
#'    Konen, O. Flasch, M. Friese, P. Koch, M. Zaefferer, B. Naujoks, M. Friese
#' @references
#' \url{http://www.springer.com/3-540-32026-1}
#' @keywords package
#' @seealso Main interface functions are \code{\link{spot}} and \code{\link{spotOptim}}.
#' Also, a graphical interface can be used with \code{\link{spotGui}}
#' @import emoa
#' @import rpart
#' @import twiddler
#' @import rgl
#' @import AlgDesign
#' @import randomForest
#' @import mco
#' @import rsm
#' @import MASS
#' 
#' @section Acknowledgments:
#' This work has been partially supported by the Federal Ministry of Education
#' and Research (BMBF) under the grants CIMO (FKZ 17002X11) and
#' MCIOP (FKZ 17N0311).
#'
#' @section Maintainer:
#' Martin Zaefferer \email{martin.zaefferer@@gmx.de} 
#End of Package Description
NA #NULL, ends description without hiding first function
###################################################################################
## spot.R consists of three parts: 
## - PART ONE: some help functions
## - PART TWO: the steps implemented as functions too
## - PART THREE: the main SPO algorithm
###################################################################################

###################################################################################
###################################################################################
# PART ONE: help functions 
###################################################################################
###################################################################################

###################################################################################
# SPOT Prepare System - loads all required packages for SPOT
# 
# installs and loads all packages that are needed for the core functionality of SPOT
# (hard coded in the function). All user defined Plugins needs to call
# \code{\link{spotInstAndLoadPackages}} to add their dependencies properly\cr
# This function is only provided for use in non-packaged version, for all packages 
# are listed in the "Depends line" of DESCRIPTION 
# 
#
# @export
# @keywords internal
###################################################################################
#spotPrepareSystem <- function(){
	### check whether necessary packages are installed and install missing packages
	# see also: depends and suggests in DESCRIPTION file.
	#necessaryPackages = c('rpart', 'emoa' ); 
	#spotInstAndLoadPackages(necessaryPackages)	
	###### default packages with various use 
	# 'emoa' - used in various functions for multi objective optimization purpose, but mainly in spotGenerateSequentialDesign
	###### default packages that are specified to be used in:  
	# spotPreditTree AND spotPlotBst: 'rpart'
	###### deleted because use was not found:  
	# 'colorspace',Color Space Manipulation 
	# 'vcd' Visualizing Categorical Data - used in ???
	# 'stats' -  R statistical functions used in ??? 
	# 'DoE.base' used in ???
	# 'car' used  in ???
	######
	## deleted from list and moved to the calling functions: 
	## rsm, tgp, randomForest, mlegp, FrF2, DoE.wrapper, AlgDesign, lhs, fields
	# spotPredictLm: 'rsm'
	# spotPredictTgp:  'tgp',
	# spotPredictRandomForest: 'randomForest',
	# spotPredictMlegp: 'mlegp',
	# spotCreateDesignFrF2 : 'FrF2',  'DoE.wrapper',
	# spotPredictDiceKriging: ,'DiceKriging' # depreciated
	# spotCreateDesignBasicDoe: 'AlgDesign'	
#}#end spotPrepareSystem

###################################################################################
## SPOT Prepare
###################################################################################
#' Prepares the configuration (spotConfig) for SPOT  
#' 
#' Set some globally important parameters for SPOT to run, creating the
#' parameter list (spotConfig) for SPOT from given usersConfigFile  (.conf -file)
#' 
#' @note For developers: this function also manages the include of all functions needed - 
#' in the packaged version this is already done when package is installed.
#'
#' @param srcPath 			the absolute path to the SPOT sources
#' @param configFile 		the absolute path including file-specifier
#' @param spotConfigUser	a list of parameters used to configure spot, usually  spotConfigUser=NA will be passed to this function, which means the configuration will only be read from the \code{configFile}, not given by manual user input. 		
#'							Notice that parameters given in spotConfigUser will overwrite both default values assigned by SPOT, AND values defined in the config file
#'							However, values not passed by spotConfigUser will still be used as defaults. If you want to see those defaults, look at \code{\link{spotGetOptions}}  
#' @return list \code{spotConfig} \cr
#' - \code{spotConfig} is the list of spot parameters created by this function
#'
#' @seealso \code{\link{SPOT}} \code{\link{spotGetOptions}} \code{\link{spot}}
#' @export
#' @keywords internal
###################################################################################
spotPrepare <- function(srcPath,configFile,spotConfigUser){
	# Close graphic windows	
	graphics.off()
	######################################
	### Load sources
	######################################	
	## Add path to files
	## everything happens relative to users configuration file
	if(file.exists(configFile)){
		setwd(dirname(configFile))
	}
	## Call configuration program that extracts infos from userconf	
	spotConfig <- spotGetOptions(srcPath=srcPath,configFile)
	## MZ 04.09.2010: New feature implemented, so user can set options in commandline when calling  spot()
	if(is.list(spotConfigUser)){
		spotConfig <- append(spotConfigUser,spotConfig)
		spotConfig <- spotConfig[!duplicated(names(spotConfig))]#Commandline Input from user will overwrite configfile/default parameters here !!
		if(file.exists(spotConfig$io.roiFileName)){ #Read in the roi, just in case that spotConfigUser contained a new roi file name
			spotConfig$alg.roi <- spotReadRoi(spotConfig$io.roiFileName,spotConfig$io.columnSep,spotConfig$io.verbosity)
			spotConfig$alg.aroi <- spotConfig$alg.roi
		}
		colnames(spotConfig$alg.roi) <- c("lower","upper","type")		
		colnames(spotConfig$alg.aroi) <- c("lower","upper","type")
	}
	if(is.function(spotConfig$alg.func)){
		spotConfig$alg.tar.func<-spotConfig$alg.func
		spotConfig$alg.func<-"spotOptimInterface"		
	}
	else if(!is.character(spotConfig$alg.func)){
		stop("The optimization target function is neither a character string, nor a function handle")
	}	
	## MZ 30.08.2012: Continue stopped or broken SPOT runs (broken runs only continued if file mode enabled)
	if(spotConfig$spot.continue){
		if(spotConfig$spot.fileMode && file.exists(spotConfig$io.resFileName)) {
			spotConfig$alg.currentResult <- spotGetRawResData(spotConfig)$rawD
		}	
	}
	
#	if (spotConfig$spot.ocba == TRUE){#Bugfix: If ocba is chosen, makes sure that max repeats, and initial repeats are more than 1. However this will still crash with noise=0
#		if (!is.na(spotConfig$init.design.repeats)){
#			if (spotConfig$init.design.repeats <= 1){
#				spotConfig$init.design.repeats=2
#			}
#		}
#		if (!is.na(spotConfig$seq.design.maxRepeats)){
#			if (spotConfig$seq.design.maxRepeats <= 1){
#				spotConfig$seq.design.maxRepeats=2
#			}
#		}
#	}
	class(spotConfig)<-"spotConfig" #TODO class might yield slow-down!
	spotConfig
} # end spotPrepare()

###################################################################################
###################################################################################
## PART TWO:  The SPO Steps 
###################################################################################
###################################################################################

###################################################################################
#' SPOT Step: Initialize (First SPOT- Step)
#' 
#' Creates a sequential design based on the results derived so far. Therefor it is
#' essential to have another design evaluated before and have a .res file to use.
#' afterwards the design is extended by 4 columns: CONFIG, REPEATS,STEP, SEED 
#'
#' uses the functions \code{spotConfig$init.design.func} and \code{link{spotWriteDes}}
#' that writes a  design to the file <xxx>.des
#' 
#' @param spotConfig the list of all parameters is given, but the used ones are: \cr 
#'   \code{spotConfig$init.design.func} holds the spotCreateDesign<XXX> function to be used 
#'			for building an initial design. \cr 
#'   \code{spotConfig$init.design.size} number of points that should be created for the initial design \cr
#'   \code{spotConfig$init.design.retries} gives the number of trials to find a design with the greatest distance between points, (default is 1)\cr 
#'   \code{spotConfig$init.design.repeats} number of repeats for one initial design-point\cr
#'   \code{spotConfig$alg.seed} seed value for reproducible runs\cr
#'   \code{spotConfig$srcPath} source path as given when spot() is called (or uses default)\cr
#'   \code{spotConfig$io.verbosity} verbosity for command window output, which is passed to the output function
#' @export
###################################################################################
spotStepInitial <- function(spotConfig) {
  	## Sets the seed for all random number generators in SPOT
	set.seed(spotConfig$spot.seed) 
	#clear old  data 
	spotConfig$alg.currentResult<-NULL
	spotConfig$alg.currentBest<-NULL
	
	spotWriteLines(spotConfig$io.verbosity,2,"Create Inital Design", con=stderr());
	if(!exists(spotConfig$init.design.func))stop(paste("The design function name", spotConfig$init.design.func, "is not found in the workspace \n
		Please make sure to load the design function in the workspace, or specify the correct function in spotConfig$init.design.func" ))
	##
	## write actual region of interest file (same data as roi file)	
	## TODO: Add type information to aroi file
	A <- spotConfig$alg.roi	
	A <- cbind(row.names(A), A)  
	colnames(A) <- c("name", "lower", "upper", "type")	
	if(spotConfig$spot.fileMode){
		spotWriteAroi(A,spotConfig$io.verbosity,spotConfig$io.columnSep,spotConfig$io.aroiFileName)	
	}
	spotConfig$alg.aroi<-spotConfig$alg.roi
	if(spotConfig$init.design.size>0){
		initDes<-eval(call(spotConfig$init.design.func, 
					spotConfig,
					spotConfig$init.design.size,
					spotConfig$init.design.retries))
	}else{
		initDes <- NULL
	}
	
	#add manually specified design points
	if(!is.null(spotConfig$init.design.man)){
		colnames(spotConfig$init.design.man) = rownames(spotConfig$alg.roi)
		initDes <- rbind(initDes,spotConfig$init.design.man)
	}
	
	if(is.null(initDes)){
		stop("Initial Design for SPOT is empty. Set spotConfig$init.design.size to a value larger than zero, or specify design points manually in spotConfig$init.design.man.")
	}
	
	## FIRST COLUMN ADDED: Named "CONFIG" - holding a count variable: 
	## number of the configuration provided
	configNumber<-1:nrow(initDes)
	initDes <- cbind(initDes,configNumber)
	colnames(initDes)[ncol(initDes)] <- "CONFIG"
	
	## SECOND COLUMN ADDED: 
	## number of repeats for the initial design points as "repeats" 
	initDes <- cbind(initDes,spotConfig$init.design.repeats)
	colnames(initDes)[ncol(initDes)] <- "REPEATS"
	## THIRD COLUMN ADDED: column documenting the number of configurations so far (steps-column)
	## initially the number of steps is 0 (refers to auto.loop.steps)
	initDes <- cbind(initDes,0)
	colnames(initDes)[ncol(initDes)] <- "STEP"

	## FORTH COLUMN ADDED: 
	## Named "SEED" - holding the number of the seed for the randomgenerator
	## used (same seed provides reproducable runs)
	seed <- spotConfig$alg.seed
	## could be considering the last used seed according to the last res, 
	## but not yet considered here	
	initDes <- cbind(initDes,seed)
	colnames(initDes)[ncol(initDes)] <- "SEED"	
	if (spotConfig$spot.fileMode){
		if (file.exists(spotConfig$io.desFileName)){
			file.remove(spotConfig$io.desFileName)
		}
		## write the design to a NEW .des-file 
		spotWriteDes(initDes,spotConfig$io.verbosity,spotConfig$io.columnSep,spotConfig$io.desFileName)		
		## Now delete the old .res and .bst files
		if (spotConfig$init.delete.previous.files & file.exists(spotConfig$io.bstFileName)){
			file.remove(spotConfig$io.bstFileName)
		}
		if (spotConfig$init.delete.previous.files & file.exists(spotConfig$io.resFileName)){
			file.remove(spotConfig$io.resFileName)
		}
	}
	spotConfig$alg.currentDesign<-initDes		
	spotConfig
}

###################################################################################
## Second Step: Algorithm Call 
###################################################################################
#' SPOT Step Algorithm Call 
#'
#' This is the second SPOT Step after step "initial" - but also needed 
#' after each step "sequential", and is a call frame for the algorithm-call.
#'
#' The algorithm is the heart of what the user must provide, but SPOT should be 
#' able to handle them in the most flexible manner. This function is an interface to the algorithm,
#' given as a R-function.
#' 
#' @param spotConfig the list of all configuration parameters, but most important ones are:\cr 
#'   \code{spotConfig$alg.func} the name of the R target function \cr
#'   \code{spotConfig$io.apdFileName} filename for the problem definition of the algorithm, 
#' 			first parameter of the generically defined R-function spotConfig$alg.func \cr
#'   \code{spotConfig$io.desFileName} filename for the input of the algorithm, 
#' 			second parameter of the generically defined R-function spotConfig$alg.func \cr
#'   \code{spotConfig$io.resFileName} filename for the output of the algorithm
#' 			third parameter of the generically defined R-function spotConfig$alg.func\cr
#'   \code{spotConfig$io.verbosity} verbosity for command window output, which is passed to the output function
#' @param ... additional parameters to be passed on to target function which is called inside alg.func
#'
#' @seealso \code{\link{SPOT}} \code{\link{spot}} \code{\link{spotStepInitial}}
#' \code{\link{spotStepSequential}}
#' @export
####################################################################################
spotStepRunAlg <- function(spotConfig,...){
	spotWriteLines(spotConfig$io.verbosity,2,paste("spotStepRunAlg started with ",spotConfig$alg.func,sep=""))
	if(!exists(spotConfig$alg.func))stop(paste("The target function name", spotConfig$alg.func, "is not found in the workspace \n
		Please make sure to load the target function in the workspace, or specify the correct function in spotConfig$alg.func"))
	#browser()
	#spotConfig<-eval(call(spotConfig$alg.func, spotConfig,...)) 
	spotConfig<-do.call(spotConfig$alg.func, args=list(spotConfig,...)) #this allows further arguments dot-dot-dot
}

###################################################################################
## Third Step: Sequential
#' SPOT Step Sequential
#'
#' Third SPOT Step to generate a sequential new design, this
#' is mainly a call of \code{\link{spotGenerateSequentialDesign}}
#' 
#' Creates a sequential design based on the results derived so far. Therefor it is
#' essential to have another design evaluated before and have a .res file to use.
#' It uses the functions \code{\link{spotGenerateSequentialDesign}} and \code{\link{spotWriteDes}}
#' writes a sequential design to the file <xxx>.des
#' 
#' @param spotConfig the list of all parameters is given, but the used ones are: \cr
#'   \code{spotConfig$io.resFileName} is checked for existence is not, function fails with error\cr
#'   \code{spotConfig$algSourceSrcPath} needed for the error message \cr
#'   \code{spotConfig$userConfFileName} needed for the error message\cr
#' @export
###################################################################################
spotStepSequential <- function(spotConfig) {
	spotWriteLines(spotConfig$io.verbosity,2,"Create Sequential Design", con=stderr())
	if(spotConfig$spot.fileMode){
		if (!file.exists(spotConfig$io.resFileName)){
			stop("Error in spot.R::spotStepSequential:
			.res file not found, spotStepAlgRun() has to be executed before.")		
		}
	}else{
		if(!nrow(spotConfig$alg.currentResult)>0){
			stop("Error in spot.R::spotStepSequential:
			result data not found, spotStepAlgRun() has to be executed before.")		
		}
	}
	
	##NOTE: the following code was moved to spotGenerateSequentialDesign, when merging with ocba
	########MZ: Now first check for var of the y-values. Only if y varies (i.e. function is noisy and evaluations are repeated) use ocba
	# varies=TRUE;
	# if(spotConfig$spot.ocba == TRUE){ # only needs to be checked for ocba=TRUE
		# if(spotConfig$spot.fileMode){
			# res<-spotGetRawResData(spotConfig)
			# spotConfig<-res$conf
			# rawData<-res$rawD
			# res<-NULL
		# }else{
		  # rawData=spotConfig$alg.currentResult;
		# } 
		# z <- split(rawData[,spotConfig$alg.resultColumn], rawData$CONFIG);
		#########varY <- sapply(as.data.frame(z),var);
		# varY <- sapply(z,var);
		# for (i in 1:length(varY)){
			# if (is.na(varY[i])||is.nan(varY[i])||is.null(varY[i])||(varY[i]==0)){
				# varies=FALSE;
			# }
		# }
	# }
	############Now call sequential design function
	# if ((spotConfig$spot.ocba == TRUE)&(varies == TRUE)){
		# spotConfig <- spotGenerateSequentialDesignOcba(spotConfig);	
    # }
    # else if (spotConfig$spot.ocba == FALSE) {
		# spotConfig <- spotGenerateSequentialDesign(spotConfig);
    # }
	# else{
		# stop("
		
		# There is no variance for some point(s) in the current design.
		
		# Therefore OCBA cannot be used. Possible reasons are a target 
		# function without noise, or the design points are not repeated.
		# SPOT with OCBA makes only sense if the target function is noisy.
		# If a non noisy function is used, the default settings should 
		# be adopted,	as described in the help of spot() or spotOptim().
		# That means: either use spot.ocba=FALSE, or set the repeats
		# (init.design.repeats) to values larger
		# than 1.
		
		# The current variance vector is for the used design points is: 
		# ",paste(varY," "))
	# }
	spotConfig <- spotGenerateSequentialDesign(spotConfig)
}

###################################################################################
## Forth Step Report
###################################################################################
#' SPOT Step Report
#'
#' Forth and last step for SPOT, that is by default a call of \link{spotReportDefault}
#' 
#' This step provides a very basic report about the .res-file, based on settings in the \code{spotConfig}
#' The mainly used parameters of \code{spotConfig} is \code{spotConfig$report.func},
#' specifying which report shall be called. The user can specify his own report and should set the 
#' value {report.func} in the configuration file according to the specification rules
#' given. If nothing is set, the default report is used.
#' 
#' @param spotConfig the list of all parameters is given, it is forwarded to the call of the report-function
#' @seealso \code{\link{SPOT}} \code{\link{spot}} \link{spotReportDefault} \code{\link{spotGetOptions}} 
#' 
#' @export
###################################################################################
spotStepReport <- function(spotConfig) {
	if(!exists(spotConfig$report.func))stop(paste("The report function name", spotConfig$report.func, "is not found in the workspace \n
		Please make sure to load the report function in the workspace, or specify the correct function in spotConfig$report.func" ))	
	if(is.null(spotConfig$alg.currentResult))spotConfig$alg.currentResult<- spotGetRawResData(spotConfig)$rawD;	
	spotConfig<-eval(call(spotConfig$report.func, spotConfig))
}

###################################################################################
## Step Auto
###################################################################################
#' SPOT Step Auto Opt
#'
#' spotStepAutoOpt is the default task called, when spot is started.  
#' 
#' The \code{auto} task calls the tasks \code{init} and \code{run} once
#' and loops \code{auto.loop.steps} times over the steps \code{seq}  and \code{run}
#' finalising the function with a call of the report function. Instead of \code{auto.loop.steps}
#' also \code{auto.loop.nevals} can be used as a stopping criterion.
#' 
#' @param spotConfig the list of all parameters is given, it is forwarded to the call of the report-function
#' the used parameters of spotConfig are just spotConfig$auto.loop.steps
#' specifying the number of meta models that should be calculated
#' @param ... additional parameters to be passed on to target function which is called inside alg.func
#' 
#' @seealso \code{\link{SPOT}} \code{\link{spot}} \code{\link{spotStepInitial}}
#' \code{\link{spotStepSequential}} \code{\link{spotStepRunAlg}} \code{\link{spotStepReport}} 
#' \code{\link{spotGetOptions}} 
#' @export
###################################################################################
spotStepAutoOpt <- function(spotConfig,...){ 
	if(!spotConfig$spot.continue || is.null(spotConfig$alg.currentResult)){
		spotConfig=spotStepInitial(spotConfig);
		spotConfig=spotStepRunAlg(spotConfig,...)
	}
	j <- max(spotConfig$alg.currentResult$STEP)
	k <- nrow(spotGetRawDataMatrixB(spotConfig));
	if(!is.null(spotConfig$spot.catch.error)){
		res<-tryCatch({ #This function will catch crashes and interrupts, but still return the last valid spotConfig list, to recover any available results and settings.
		while (j <= spotConfig$auto.loop.steps && k < spotConfig$auto.loop.nevals){	
			spotWriteLines(spotConfig$io.verbosity,2,paste("SPOT Step:", j), con=stderr());
			spotConfig=spotStepSequential(spotConfig);		
			spotConfig=spotStepRunAlg(spotConfig,...);
			k <- nrow(spotGetRawDataMatrixB(spotConfig));	
			j <- j+1;}
			}, interrupt = function(ex) {			
				cat("An interrupt was detected in spotStepAutoOpt.\n");
				print(ex);
			}, 	error = function(ex) {
				cat("An error was detected in spotStepAutoOpt.\n");
				print(ex);				
			}, finally = {}
		) #tryCatch end.	
	}else{	
		res<-tryCatch(
			{ #This function will catch crashes and interrupts, but still return the last valid spotConfig list, to recover any available results and settings.
				while (j <= spotConfig$auto.loop.steps && k < spotConfig$auto.loop.nevals){	
					spotWriteLines(spotConfig$io.verbosity,2,paste("SPOT Step:", j), con=stderr())
					spotConfig=spotStepSequential(spotConfig)
					spotConfig=spotStepRunAlg(spotConfig,...)
					k <- nrow(spotGetRawDataMatrixB(spotConfig))
					j <- j+1
				}
			}, interrupt = function(ex) {			
				cat("An interrupt was detected in spotStepAutoOpt.\n")
				print(ex)
			}, finally = {}
		) #tryCatch end.
	}	
	if(!is.null(res)){#&&  any(class(res)=="interrupt")){
		cat("A crash or interrupt ocurred, most recent spotConfig list is returned, to allow recovery of results.\n")
		return(spotConfig)
	}
	if(spotConfig$io.verbosity>2){  
		mergedData <- spotPrepareData(spotConfig)
		spotConfig=spotWriteBest(mergedData, spotConfig)
		spotPlotBst(spotConfig)   
	} 
	spotConfig=spotStepReport(spotConfig)
	spotConfig
}
###################################################################################
## Step Meta
###################################################################################
#' SPOT Step Meta 
#'
#' Attention: This feature is work in progress, documentation is not up to date.
#'
#' The \code{meta} task calls spotStepMetaOpt which itself calls \code{\link{spot}}
#' with several different fixed
#' parameters to provide a mixed optimization mechanism: analyse a fully qualified 
#' test of some parameters and the intelligent optimization of other parameters.
#' e.g. the number of the dimension of a problem etc.
#' 
#' To start this step you could for example do this:\cr
#' \code{spot("configFileName.conf","meta")}\cr 
#' 
#' @param spotConfig the list of all parameters is given
#' 
#' @seealso \code{\link{spotGetOptions}} 
#' @export
###################################################################################
spotStepMetaOpt <- function(spotConfig) {
	#spotInstAndLoadPackages("AlgDesign")
	if(is.null(spotConfig$report.meta.func))spotConfig$report.meta.func = "spotReportMetaDefault";	
	### Delete old FBS file
	if(file.exists(spotConfig$io.fbsFileName)) {
		unlink(spotConfig$io.fbsFileName)
	} 
		
	myList<-spotConfig$meta.list
	mySetList<-spotConfig$meta.conf.list
	nVars<-length(myList)
	x <- as.numeric(lapply(myList, length))	# create a vector "x" holding the length of each variable
	if (nVars==1){	# full factorial design with indicies for all combinations:
		dat <- matrix(1:x, byrow = TRUE)
	}
	else{
		dat<- gen.factorial(x,varNames=names(myList),factors="all")
	}	
	for (j in 1:nrow(dat)) {## Loop over full factorial combinations of all parameters specified in .meta
		graphics.off() ## close all remaining graphic devices - from old spotStepAutoOpt Runs
		myFbs<-list()
		newConf<-list()
		newConfSet<-list()
		for (k in 1:nVars) {
			# left side of the  assignment 
			## the factorial value of the kth variable for this dat[j]-row is assigned to a character variable:
			newConf[[names(myList[k])]] <-  myList[[k]][[dat[j,k]]]
			for (ii in 1:length(mySetList[[k]])){
				if(length(mySetList[[k]])>0){
					newConfSet[[names(mySetList[[k]][ii])]]<-mySetList[[k]][[ii]][dat[j,k]]
				}
			}
		}
		myFbs <- c(myFbs, newConf)
		newConf<-append(newConf,newConfSet)
		newSpotConfig<-append(newConf,spotConfig)	## create a temporary spotConfig for the calling of spotStepAuto 
		newSpotConfig$spot.fileMode=FALSE;
		newSpotConfig <- newSpotConfig[!duplicated(names(newSpotConfig))]; ## delete unneeded entries
		newSpotConfig=spot(spotConfig=newSpotConfig) ##  THIS calls spot for ONE configuration of the meta run
		tmpBst<-newSpotConfig$alg.currentBest;
		design = as.list(dat[j,])
		names(design)=paste(names(design),"NUM", sep="")
		myFbsFlattened <- spotMetaFlattenFbsRow(append(myFbs,design))
	
		dataTHIS<-as.data.frame(cbind(tmpBst[nrow(tmpBst),],myFbsFlattened))
		
		if(file.exists(spotConfig$io.fbsFileName)) {
			dataLAST<-as.data.frame(read.table(file=spotConfig$io.fbsFileName,header=TRUE))
			data<-merge(dataLAST,dataTHIS,all=TRUE,sort=FALSE)
		}
		else{
			data=dataTHIS
		}
		write.table(file=spotConfig$io.fbsFileName,
				data,
				row.names = FALSE,
				col.names = TRUE,
				sep = " ",
				append = FALSE,
				quote=FALSE)
				
	} # for (j in 1:nrow(dat))... (loop over full factorial design)
	spotConfig$meta.fbs.result=data
	if(!exists(spotConfig$report.meta.func))stop(paste("The meta report function name", spotConfig$report.meta.func, "is not found in the workspace \n
		Please make sure to load the meta report function in the workspace, or specify the correct function in spotConfig$report.meta.func"))	
	spotConfig<-eval(call(spotConfig$report.meta.func, spotConfig))
}

############# end function definitions ############################################################

###################################################################################################
## PART THREE: SPOT: The Program
###################################################################################################
#' Main function for the use of SPOT
#' 
#' Sequential Parameter Optimization Toolbox (SPOT) provides a toolbox for the 
#' sequential optimization of parameter driven tasks. 
#' Use \code{\link{spotOptim}} for a \code{\link{optim}} like interface
#'
#' The path given with the \code{userConfigFile} also fixes the working directory used
#' throughout the run of all SPOT functions. All files that are needed for input/output
#' can and will be given relative to the path of the userConfigFile (this also holds for 
#' the binary of the algorithm). This refers to files that are specified in the configFile
#' by the user. 
#'
#' It is of major importance to understand that spot by default expects to optimize noisy functions. That means, the default settings of spot,
#' which are also used in spotOptim, include repeats of the initial and sequentially created design points. Also, as a default OCBA
#' is used to spread the design points for optimal usage of the function evaluation budget. OCBA will not work when there is no variance in the data.
#' So if the user wants to optimize non-noisy functions, the following settings should be used:\cr
#' \code{spotConfig$spot.ocba <- FALSE}\cr
#' \code{spotConfig$seq.design.maxRepeats <- 1}\cr
#' \code{spotConfig$init.design.repeats <- 1}\cr
#'
#' @param configFile	the absolute path including file-specifier, there is no default, this value should always be given
#' @param spotTask		[init|seq|run|auto|rep] the switch for the tool used, default is "auto" 
#' @param srcPath		the absolute path to user written sources that extend SPOT, the default(NA) will search for sources in the path <.libPath()>/SPOT/R  
#' @param spotConfig	a list of parameters used to configure spot, default is spotConfig=NA, which means the configuration will only be read from the \code{configFile}, not given by manual user input. 		
#'						Notice that parameters given in spotConfig will overwrite both default values assigned by SPOT, AND values defined in the Config file
#'						However, values not passed by spotConfig will still be used as defaults. If you want to see those defaults, look at \code{\link{spotGetOptions}}  
#' @param ... additional parameters to be passed on to target function which is called inside alg.func. Only relevant for spotTask "auto" and "run".
#' @note \code{spot()} expects char vectors as input, e.g. \code{spot("c:/configfile.conf","auto")}
#' @seealso \code{\link{SPOT}}, \code{\link{spotOptim}}, \code{\link{spotStepAutoOpt}},  \code{\link{spotStepInitial}},
#' \code{\link{spotStepSequential}}, \code{\link{spotStepRunAlg}}, \code{\link{spotStepReport}} 
#' @export
###################################################################################################
spot <- function(configFile="NULL",spotTask="auto",srcPath=NA,spotConfig=NA,...){
	writeLines("spot.R::spot started ") #bugfix MZ: spotWriteLines will not allways work here, since spotConfig could be NA
	callingDirectory<-getwd()
	if(!(configFile=="NULL")&!file.exists(configFile)){
		stop("Error, configFile not found (or not \"NULL\")")
	}	
	if(is.na(srcPath)){
		for(k in 1:length(.libPaths())){ 
			if(file.exists(paste(.libPaths()[k],"SPOT","R",sep="/"))){
				srcPath<-(paste(.libPaths()[k],"SPOT","R",sep="/"))
				break;
			}
		}
	}
	## PRELIMINARIES 1: load all functions belonging to SPOT - not necessary if provided SPOT is installed as package - useful for developers...
	spotConfig<-spotPrepare(srcPath,configFile,spotConfig)

	## SWITCH task according to the extracted from command line 
	resSwitch <- switch(spotTask
			, init=, initial=spotStepInitial(spotConfig) # First Step
			, seq=, sequential=spotStepSequential(spotConfig) # Second Step
			, run=, runalg=spotStepRunAlg(spotConfig,...)	# Third Step
			, rep=, report=spotStepReport(spotConfig)		# Fourth Step
			, auto=, automatic=spotStepAutoOpt(spotConfig,...)	# Automatically call First to Forth Step
			, meta=spotStepMetaOpt(spotConfig)	# Automatically call several spotStepAutoOpt - Runs to provide a systematic testing tool an fixed Parameters in .apd file	
			, "invalid switch" # return this at wrong CMD task
	);
	## ERROR handling 
	## valid switch returns null, otherwise show error warning and short help
	if (is.character(resSwitch) && resSwitch == "invalid switch") {
		#spotWriteLines(spotConfig$io.verbosity,0,paste("ERROR, unknown task:", spotTask), con=stderr());
		#spotWriteLines(spotConfig$io.verbosity,0,"\nValid tasks are:\
		#				auto       - run tuning in automated mode\
		#				initial    - to create an initial design\
		#				run        - start the program, algorithm, simulator\
		#				sequential - to create further design points\
		#				report     - to generate a report from your results"
		#		, con=stderr());
		stop(paste("ERROR, unknown task:", spotTask, "\nValid tasks are:\
						auto       - run tuning in automated mode\
						initial    - to create an initial design\
						run        - start the program, algorithm, simulator\
						sequential - to create further design points\
						report     - to generate a report from your results" ))		
	}
	# go back to -  well where ever you came from
	setwd(callingDirectory)
	resSwitch
}

###################################################################################################
#' Print function for spotConfig class 
#'
#' Print function to summarize a spotConfig.
#'
#' @rdname print
#' @method print spotConfig
# @S3method print spotConfig
#' @param x	spotConfig
#' @param ... additional parameters	
#' @export
#' @keywords internal
#####################################################################################
print.spotConfig <- function(x, ...){ 	#Remark: this overwrites the print method for the spotConfig class
										#This class is set for spotConfig only once, in spotPrepare.
	writeLines(paste("This is a spotConfig list"))
	writeLines("Current list content:")
	print(names(x))
	writeLines("Use \"listname[]\" to print all list values, e.g. spotConfig[].")
	writeLines("See the help of spotGetOptions for more information.")
} 
