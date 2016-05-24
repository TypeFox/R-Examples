## Experimental research in evolutionary computation
## author: thomas.bartz-beielstein@fh-koeln.de
## http://www.springer.com/3-540-32026-1
##
## Copyright (C) 2009 T. Bartz-Beielstein and C. Lasarczyk
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
###################################################################################################
#' Set all options by conf-file or by default
#' 
#' spotGetOptions \cr
#' 1.) sets default values \cr 
#' 2.) overwrites all default values by the settings the user provides with the config file (.conf-file) \cr
#' All options described here, that are not marked as "internal variable" may be changed 
#' by the user. This will be done by reading the ".conf"-file that the user has specified 
#' as the first (and maybe sole) parameter to the function spot(). 
#' To change this default value of a variable, simply write a line into the ".conf"-file following
#' this syntax:\cr
#' <variable>=<value> e.g.: \code{spot.seed=54321}\cr
#' This function will do even more: the user may define his own variables in the .conf-file and may use them
#' in user written plugins. All plugins will get the whole list of options with the parameter 
#' "spotConfig". As a result a variable given in the .conf file as \cr 
#' \code{my.var=37} \cr
#' may be referred to by spotConfig$my.var and can be used in all functions - especially in the functions 
#' that are designed to be open to adaptions where ever necessary. 
#' 
#' @param srcPath the absolute path to the SPOT sources
#' @param configFileName users config file (.conf) the absolute path including file-specifier of the user config File
#' @return spotGetOptions returns the list of all SPOT options  created by this function:
#' 			\item{auto.loop.steps}{[\code{Inf}] number of iterations the loop over all SPOT-steps should be repeated}
#' 			\item{auto.loop.nevals}{[\code{200}] budget  of algorithm/simulator runs
#' 			 - most important parameter for run-time of the algorithm in case the spot-function is called with the "auto"-task }
#' 			\item{spot.continue}{[\code{FALSE}] boolean, SPOT will try to continue based on existing results in spotConfig or .res file if this value is TRUE}
#' 			\item{spot.fileMode}{[\code{TRUE}] boolean, that defines if files are used to read and write results (which is the "classic" spot procedure) or if SPOT will only use the workspace to store variables.}
#' 			\item{spot.seed}{[\code{123}] global seed setting for all random generator dependent calls within SPOT. same seed shall repeat same results, \cr
#' 						BUT: please note: this is NOT the seed for the algorithm! see alg.seed}
#' 			\item{alg.func}{[\code{"spotFuncStartBranin"}] target function to be optimized by SPOT.\cr SPOT searches for the given function name in workspace}
#' 			\item{}{If alg.func is a string, SPOT expects an interface like \code{\link{spotFuncStartBranin}}.\cr
#'					If alg.func is a function, SPOT expects a function of type y=f(x)\cr
#' 					(see: \code{\link{spotOptimInterface}}).} 
#' 			\item{alg.resultColumn}{[\code{"Y"}] string to indicate the name of the result column. This might be automatically reset if \code{\link{spotOptimInterface}} is used, i.e. if a function is passed to alg.func}
#' 			\item{}{It can be a vector of strings, if multi objective optimization is done.}
#'			\item{alg.seed}{[\code{1234}] seed for random generator to be used by the user defined algorithm. 
#'						This is needed to reproduce the results. Set to NA if seed should not be reset. }
#' 			\item{alg.roi}{internal parameter for the initial region of interest (do not try to set this one, it will be overwritten with default values).
#'                       It is used to provide an easy to use matrix with the data from the  ".roi"-file (= Region Of Interest)}
#'			\item{}{The user can create an ROI matrix with the \code{\link{spotROI}} constructor.}
#' 			\item{alg.aroi}{internal parameter for the actual region of interest (do not try to set this one, it will be overwritten with default values). 
#' 							It is used to provide an easy to use matrix with the data from the  ".aroi"-file (= Actual Region Of Interest)}
#' 			\item{alg.currentDesign}{[usually not changed by user] data frame of the design that will be evaluated by the next call to \link{spotStepRunAlg}} 
#' 			\item{alg.currentResult}{[usually not changed by user] data frame that contains the results of the target algorithm runs} 
#' 			\item{alg.currentBest}{[usually not changed by user] data frame that contains the best results of each step conducted by spot} 
#' 			\item{io.columnSep}{[\code{""}] column separator for the input/output files, default means: arbitrary whitespace sequence, 
#' 							should be set by the value you want to have between your columns}
#' 			\item{io.apdFileName}{[depends: \code{<configFileName>.apd}] name of the .apd -file (Algorithm Problem Definition file, holding all specification the user written algorithm needs to perform a complete optimization)}
#' 			\item{io.roiFileName}{[depends: \code{<configFileName>.roi}] name of the .roi -file (Region Of Interest - File, holding all varying parameters and constraints)}
#' 			\item{io.desFileName}{[depends: \code{<configFileName>.des}] name of the .des -file (DESign file, the file the user written algorithm uses as input to the parameters it should change)}
#' 			\item{io.resFileName}{[depends: \code{<configFileName>.res}] name of the .res -file (RESult file) the user written algorithm has to write its results into this file }
#' 			\item{io.bstFileName}{[depends: \code{<configFileName>.bst}] name of the .bst -file (BeST file) the result-file will be condensed to this file }
#' 			\item{io.pdfFileName}{[depends: \code{<configFileName>.pdf}] name of the .pdf -file the default report will write its summary of results in this pdf file }
#' 			\item{io.fbsFileName}{[depends: \code{<configFileName>.bst}] name of the .fbs -file (Final BestSolution file) collects all final best values of all .bst files during a .meta-run }
#' 			\item{io.verbosity}{[\code{3}] level of verbosity of the program, 0 should be silent and 3 should produce all output- sometimes just interesting for the developer...}
#' 			\item{init.design.func}{[\code{"spotCreateDesignLhd"}] name of the function to create an initial design. Please also see the notes SPOT - extensions}
#' 			\item{init.design.size}{[\code{10}] number of initial design points to be created. Required by some space filling design generators. Will be used in the <init.design.func>.R-file. If \code{=NA} a value is calculated by formula.}
#' 			\item{init.design.retries}{[\code{100}] number of retries the initial designs should be retried to find randomly a design with maximum distance between the points 
#' 								This parameter will be ignored if the function is deterministic (like doe)}
#' 			\item{init.design.repeats}{[\code{2}] number of repeats for each design point to be called with the <alg.func>}
#' 			\item{init.delete.previous.files}{[\code{TRUE}] delete an existing res or bst file.  Should be set to FALSE if a SPOT run is resumed, after crash or user triggered stop.}
#' 			\item{seq.design.size}{[\code{10000}] number of sequential design points to be created}
#' 			\item{seq.design.retries}{[\code{10}] number of retries the initial designs should be retried to find randomly a design with maximum distance between the points,  
#' 								This parameter will be ignored if the function is deterministic (like doe)}
#' 			\item{seq.design.oldBest.size}{[\code{1}] number of the best already evaluated design points that should be taken into consideration for the next sequential designs (e.g. for to have an appropriate number of repeats}
#' 			\item{seq.design.new.size}{[\code{2}] according to the predictor the new design points during the seq step are ordered by their expected values. This parameter states how many new design points should be evaluated}
#' 			\item{seq.design.maxRepeats}{[\code{NA}] each design point is to be evaluated several times for statistically sound results. The number of "repeats" will increase, but will not exceed this seq.design.maxRepeats - value }
#' 			\item{seq.design.increase.func}{[\code{"spotSeqDesignIncreasePlusOne"}] functional description of how the repeats are increased (until the seq.design.maxRepeats are reached). Default increases the number of repeats by adding one.}
#' 			\item{seq.design.func}{[\code{"spotCreateDesignLhd"}] name of the function to create sequential design. Please also see the notes SPOT - extensions}
#' 			\item{seq.mco.selection}{[\code{"hypervol"}] selection scheme for new design candidates in case of multi objective optimization. "hypervol" considers contribution of each point, "tournament2" is a tournament selection. "tournament1" is not yet recommended for use.}
#' 			\item{seq.predictionModel.func}{[\code{"spotPredictRandomForest"}] name of the function calling a predictor. Default uses a Random Forest.}
#'			\item{seq.predictionOpt.func}{[\code{NA} If not NA this string will be interpreted as a function name. The function is expected to add a new setting to the sequential design. See \code{\link{spotPredictOptMulti}}}
#' 			\item{seq.merge.func}{ [\code{mean}] defines the function that merges the results from the different repeat-runs for a design. Default is to calculate the mean value.}
#' 			\item{seq.transformation.func}{[\code{I}] function for transformation of "Y" before new model is created, default: Identity function}
#' 			\item{seq.useAdaptiveRoi}{[\code{FALSE}] use region of interest adaptation}
#'			\item{report.func}{[\code{"spotReportDefault"}] name of the function providing the report ("spotReportSens","spotReport3d","spotReportContour") }
#' 			\item{report.io.screen}{[\code{FALSE}] report graphics will be printed to screen (FALSE=no, TRUE=yes)}
#' 			\item{report.io.pdf}{[\code{TRUE}] report graphics will be printed to pdf (FALSE=no, TRUE=yes)}
#' @note Additional settings can and will be written to the spotConfig list by other optional functions. See the documentation of these to functions for details.
#' @seealso  \code{SPOT} \code{spotPrepare}
#' @export
## hint for programmers: 
## ALL Variables defined in this function AND by the sourced .conf-file will be populated to the resulting list
## that is returned and used throughout the program as "spotConfig". This function just sets some defaults
## if a really temporary variable is needed for finally calculate a variable worth to be added to the list
## you must declare this "local" variable by a leading dot: e.g.  .dataPath 
##
spotGetOptions <- function( srcPath=".",configFileName) {	
	#######################################
	### Begin: Algorithm design related ###
	#######################################	
	## Specify name of the function
	## this may be 
	## 1) a function in the workspace, specified with a string giving the name of the function
	## 2) or else alg.func is directly an R object of type function (which takes x and returns y both as vectors or numerics)	
	alg.func <- "spotFuncStartBranin"	
	## Column name containing results
	alg.resultColumn <- "Y" 
	alg.seed <- 1234
	## ##########################
	## ##### Init related   #####
	## ##########################	
	## create design for initial step
	init.design.func <- "spotCreateDesignLhd"	
	## Initial number of design points, 
	## if NA, formula will be used
	## number of initial design points
	init.design.size <- 10
	## repeats for improving min-max designs (not used by doe)
	init.design.retries <- 100
	## number of repeated runs of each configuration:
	init.design.repeats <- 2
	## keep or delete existing resultfile and bstfile? (default: TRUE)
	init.delete.previous.files <- TRUE  
	## #################################
	## ##### Sequential Step related ###
	## #################################
	## Function for a summary of the results from the algorithm at each design point
	seq.merge.func <- mean	
	## function for transformation of "Y" before new model is created, default: Identitity function
	seq.transformation.func <- I   #TODO TBB doku
	## design create for sequentiel step: (must guarantee a higher search space)
	seq.design.func <- "spotCreateDesignLhd"	
	seq.design.size <- 200 #todo (auch bei init.design) vllt warnings wenns komplett daneben liegt
	seq.design.retries <- 10 #todo TBB 
	seq.design.maxRepeats <- NA 
	seq.design.increase.func <- "spotSeqDesignIncreasePlusOne" 
	## how many old (best) points shall be repeated 
	seq.design.oldBest.size <- 3  
	## how many new points shall be added 
	seq.design.new.size <- 3
	## ###################################
	## ##### Prediction Modell related ###
	## ###################################	
	seq.predictionModel.func <- "spotPredictRandomForest"	
 	seq.predictionOpt.func<-NA
	seq.predictionOpt.budget<-100
	seq.useAdaptiveRoi <- FALSE 
    seq.ocba.budget <- 3      
	seq.model.variance <- FALSE #whether model should predict variance (may be needed if seq.model.transform is not NA)
	seq.infill <- NA #expected improvement function to be used. NA means no EI #for transformations of model information, like EI, SEI, SMS-INFILL, etc
	## ###################################
	## ##### For MCO ###
	## ###################################	
	## type of pareto optimizaion selection scheme
	seq.mco.selection <- "hypervol" #"hypervol" sort by hypervolume contribution, considering known points; "tournament" tournament selection
	
	## #####################################
	## ##### Globally needed           #####
	## #####################################
	spot.seed <- 123
	spot.ocba <- TRUE 
	spot.continue <- FALSE
	## How many spot iterations should be performed?
	auto.loop.steps <- Inf
	auto.loop.nevals <- 200
	## #####################################
	## ##### Step      meta            #####
	## #####################################
	###	meta.keepAllFiles <- FALSE
	## ####################################################################################################
	## ##### IO related (files for input and output, and variables to specify formatting of IO-files  #####
	## ####################################################################################################
	## Name of the column, where the actual auto.loop.step will be stored in the res file
	## io.verbosity 0 means be quite, 3 means tell me everything
	io.verbosity<-0	#plot ausgabe nur mit interactive, bzw getrennt von verbosity
	## Separation of columns
	io.columnSep <- ""  
	## ###########################
	## ##### Report related #####
	## ###########################
	report.func <- "spotReportDefault"
	## Should graphical output be generated in a pdf File? FALSE  = NO
	report.io.pdf<-FALSE 
	## Should graphical output be generated on screen? FALSE  = NO
	report.io.screen<-TRUE 
	## reduce the number of relevant variables for the best-plot (that is also shown continuously)
	## New variable, write files or not
	spot.fileMode<-TRUE
	##########################################################################
	if(io.verbosity>0){
		writeLines("spotGetOptions... started", con=stderr())
	}
	.dataPath <- dirname(configFileName)
	if(io.verbosity>0){
        writeLines(paste("  Data Path (all experiment data are relevant to this location): ", .dataPath,collapse=""),con=stderr())
	}
	.genericFileNamePrefix <-  unlist(strsplit(basename(configFileName), ".", fixed = TRUE))[1]
	if(io.verbosity>0){
		writeLines(paste("  File name prefix: ", .genericFileNamePrefix,collapse=""),con=stderr())
	}
	io.resFileName <- paste(.genericFileNamePrefix,"res",sep=".")
	io.desFileName <- paste(.genericFileNamePrefix,"des",sep=".")
	io.bstFileName <- paste(.genericFileNamePrefix,"bst",sep=".")
	io.pdfFileName <- paste(.genericFileNamePrefix,"pdf",sep=".")	
	io.roiFileName <- paste(.genericFileNamePrefix,"roi",sep=".")
	io.aroiFileName <- paste(.genericFileNamePrefix,"aroi",sep=".")
	io.apdFileName <- paste(.genericFileNamePrefix,"apd",sep=".")
	io.metaFileName <- paste(.genericFileNamePrefix,"meta",sep=".")
	io.fbsFileName <- paste(.genericFileNamePrefix,"fbs",sep=".")	
	############################################################################
	### load user settings, this overwrite the defaults that are set up to this line
	### the basename of the configFile was used in the main 
	############################################################################
	userConfFileName  <-  basename(configFileName)
	.lsBeforeSource<-ls()
	#################################################################
	### load configuration
	#################################################################
	if(file.exists(userConfFileName)){
		source(userConfFileName, local=TRUE) # ! otherwise default values will not be overwritten
	}
	#################################################################
	### load configuration done!!!!!!!
	#################################################################
	.lsAfterSource<-ls()
	## give some warnings if NEW variables are created by the conf-file
	if (length(.lsDiff<-setdiff(.lsAfterSource,.lsBeforeSource))){
		spotWriteLines(io.verbosity,1,paste("Note: a new variable defined by conf-file (",userConfFileName,"):",.lsDiff))
	}else{
		spotWriteLines(io.verbosity,1,paste("  User conf loaded from: ",userConfFileName,collapse=""),con=stderr())
	}	
	spotWriteLines(io.verbosity,1,paste("  ResultFile Name : ", io.resFileName, collapse=""),con=stderr())
	spotWriteLines(io.verbosity,1,paste("  DesignFile Name : ", io.desFileName, collapse=""),con=stderr())
	spotWriteLines(io.verbosity,1,paste("  BestFile Name : ", io.bstFileName, collapse=""),con=stderr())
	spotWriteLines(io.verbosity,1,paste("  FinalBestSolution FbsFile Name : ", io.fbsFileName, collapse=""), con=stderr())
	## TBB: Added 26 Feb 2009:
	spotWriteLines(io.verbosity,1,paste("  pdfFile Name : ", io.pdfFileName, collapse=""),con=stderr())		
	spotWriteLines(io.verbosity,1,paste("  Load algorithm design (ROI): ", io.roiFileName, collapse=""),con=stderr())
	## alg.roi is a table that holds the data of the .roi-file for easy and quick use in some functions 
	if(file.exists(io.roiFileName)&&!(userConfFileName=="NULL")){   
		alg.roi <- spotReadRoi(io.roiFileName,io.columnSep,io.verbosity)
	}
	else if (userConfFileName=="NULL"){   #a default roi is used, will be overwritten by any user input
		alg.roi=spotROI(c(-1,-1),c(1,1)) 
		spot.fileMode=FALSE
	}
	else
	{
		alg.roi=NA
	}
	## at startup, the actual roi (alg.aroi) is the same as the initial roi (alg.roi)
	alg.aroi <- alg.roi
	spotWriteLines(io.verbosity,1,"spotGetOptions finished", con=stderr())
	## generate a list of ALL the defined variables (default AND user written = sourced by .conf-file!)
	## ls returns a list of all variables in this environment, that is: all the "local" variables (except the ones with leading dot
	.x<-ls()
	## now use sapply to generate the list of name/value of all variables. 
	sapply(.x, function (.x) { get(.x)}, USE.NAMES=TRUE)
}

