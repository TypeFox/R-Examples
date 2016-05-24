######################################################################################
# tdmDefaultsFill
#
#'   Default values for list \code{tdm}.  This list controls the tuning and unbiased evaluation phase.
#'
#'   When called with \code{tdm = tdmDefaultsFill()}, a new list \code{tdm} is created and returned.
#'   When called with \code{tdm = tdmDefaultsFill(mainFile="my.r"}, a new list \code{tdm} is created
#'   and returned, with the element mainFile set to the specified value.
#'   When called with \code{tdm = tdmDefaultsFill(tdm)}, an existing list \code{tdm} is filled with further default values.
#'   
#'   If \code{tdm$mainFunc} is missing, but \code{tdm$mainFile} exists, then \code{tdmDefaultsFill}
#'   will set \preformatted{tdm$mainFunc=sub(".r","",basename(tdm$mainFile),fixed=TRUE)}
#'
#'   @param tdm         (optional)
#'   @param mainFile    (optional) if given, create or overwrite tdm$mainFile with this value
#   @param mainCommand (deprecated) if given, create or overwrite tdm$mainCommand with this value
#'   @return \code{tdm}     the new / extended list,  where additional elements, if they are not yet def'd,  are set as: 
#'      \item{mainFile}{[NULL] if not NULL, source this file from the current dir. It should contain the definition of tdm$mainFunc.  }
#'      \item{mainFunc}{\code{sub(".r","",basename(tdm$mainFile),fixed=TRUE)}, if tdm$mainFile is set and tdm$mainFunc is NULL, else \code{"mainFunc"} 
#'              This is the name of the function called in \code{\link{tdmStartSpot}} and \code{\link{unbiasedRun}}   }
#     \item{mainCommand}{ deprecated, will be automatically set to \code{"result <- tdm$mainFunc(opts,dset=dset)}"}   }
#'      \item{unbiasedFunc}{["unbiasedRun"] which function to call for unbiased evaluation}
#'      \item{tuneMethod}{["spot"] other choices: "cmaes", "bfgs", ..., see \code{\link{tdmDispatchTuner}} }
#'      \item{nExperim}{[1]}
#'      \item{umode}{["RSUB"], one out of [ "RSUB" | "CV" | "TST" | "SP_T" ], see \code{\link{unbiasedRun}}}
#'      \item{timeMode}{[1] 1: proc time, 2: system time, 3: elapsed time (columns \code{Time.TST} and \code{Time.TRN} in \code{envT$theFinals}}
#'      \item{fileMode}{[FALSE] see "Note" section in \code{\link{tdmBigLoop}}   }
#'      \item{finalFile}{[NULL]  filename where to save \code{envT$theFinals}, only relevant for \code{tdm$fileMode==TRUE}}
#'      \item{experFile}{[NULL] filename where to append \code{envT$theFinals}, only relevant for \code{tdm$fileMode==TRUE}  }
#'      \item{filenameEnvT}{[NULL] filename where \code{\link{tdmBigLoop}} will save a small version of environment \code{envT}. If NULL, 
#'                save \code{envT} to \code{sub(".conf",".RData",tdm$runList[1])}. This RData file is written irrespective of \code{fileMode}'s value, 
#'                but only in case \code{spotStep=="auto"}. }
#'      \item{theSpotPath}{[NA] use SPOT's package version}
#'      \item{parallelCPUs}{[1] 1: sequential, >1: parallel execution with this many CPUs (package parallel)  }
#'      \item{parallelFuncs}{[NULL] in case tdm$parallelCPUs>1: a string vector with functions which are clusterExport'ed in addition
#'                to tdm$mainFunc.  }
#'      \item{path}{[NULL] where to search .conf and .apd file. If it is NULL, path is set 
#'                to the actual working directory at the time when tdmEnvTMakeNew is executed  }
#'      \item{runList}{ [NULL] a list of .conf filenames}
#'      \item{spotList}{[NULL] a list of .conf filenames}
#'      \item{stratified}{[NULL] see \code{\link{tdmSplitTestData}}  }
#'      \item{tdmPath}{[NULL] from where to source the R sources. If NULL load library TDMR instead.  }
#'      \item{test2.string}{["default cutoff"] }
#'      \item{optsVerbosity}{[0] the verbosity for the unbiased runs}
#'      \item{withParams}{[TRUE] list the columns with tuned parameter in final results  }
#'      \item{nrun}{[5] number of runs for unbiased run}
#'      \item{U.saveModel}{[FALSE] if TRUE, save the last model, which is trained in unbiasedRun, onto \code{filenameEnvT}}
#'      \item{tstCol}{["TST"] opts$TST.COL for unbiased runs (only for umode="TST") }
#'      \item{nfold}{[10] number of CV-folds for unbiased runs (only for umode="CV") }
#'      \item{TST.trnFrac}{[NULL] train set fraction (of all train-vali data),OVERWRITES opts$TST.trnFrac if not NULL. }
#'      \item{TST.valiFrac}{[NULL] validation set fraction (of all train-vali data), OVERWRITES to opts$TST.valiFrac if not NULL. }
#'      \item{TST.testFrac}{[0.2] test set fraction (of *all* data) for unbiased runs (only for umode="RSUB" or ="SP_T") }
#'      \item{CMA.propertyFile}{[NULL] (only for CMA-ES Java tuner) see \code{\link{cma_jTuner}}. }
#'      \item{CMA.populationSize}{[NULL] (only for CMA-ES Java tuner) see \code{\link{cma_jTuner}}. }
#'
#' @note 
#'      The settings \code{tdm$TST.trnFrac} and \code{tdm$TST.valiFrac} allow to set programmatically certain values for
#'      \code{opts$TST.trnFrac} and \code{opts$TST.valiFrac} *after* \code{opts} has been read from APD file. So use 
#'      \code{tdm$TST.trnFrac} and \code{tdm$TST.valiFrac} with CAUTION!
#'
#'      For \code{tdm$timeMode}, the 'user time' is the CPU time charged for the execution of user instructions of the calling process. 
#'      The 'system time' is the CPU time charged for execution by the system on behalf of the calling process. 
#'      The 'elapsed time' is the 'real' (wall-clock) time since the process was started.
#'
#' @author Wolfgang Konen, Patrick Koch, 2011 - 2013
#' @export
######################################################################################
tdmDefaultsFill <- function(tdm=NULL,mainFile=NULL) {
  if (is.null(tdm)) tdm <- list();

  if (!is.null(mainFile)) tdm$mainFile=mainFile;
  #if (!is.null(mainCommand)) tdm$mainCommand=mainCommand;
  
  if (is.null(tdm$mainFunc))  {
    tdm$mainFunc <- ifelse(is.null(tdm$mainFile), "main_sonar", sub(".R","",sub(".r","",basename(tdm$mainFile),fixed=TRUE)));   
    # e.g. tdm$mainFunc="main_myTask" if tdm$mainFile="C:/myDir/main_myTask.r"
  }
  #if (is.null(tdm$mainCommand)) tdm$mainCommand <- paste("result <- ", tdm$mainFunc,"(opts,dset=dset)",sep=" ");
  if (is.null(tdm$unbiasedFunc)) tdm$unbiasedFunc <- "unbiasedRun";
  if (is.null(tdm$tuneMethod)) tdm$tuneMethod <- "spot";
  if (is.null(tdm$nExperim)) tdm$nExperim <- 1;
  if (is.null(tdm$umode)) tdm$umode <- "RSUB";
  if (is.null(tdm$timeMode)) tdm$timeMode <- 1;             # user time
  if (is.null(tdm$fileMode)) tdm$fileMode <- FALSE;         # 10/2012 /WK/ changed from prior 'TRUE'
  if (is.null(tdm$optsVerbosity)) tdm$optsVerbosity <- 0;   # the verbosity for the unbiased runs
  if (is.null(tdm$withParams)) tdm$withParams <- TRUE;      # list the columns with tuned parameter in final results 
  if (is.null(tdm$theSpotPath)) tdm$theSpotPath <- NA;
  if (is.null(tdm$parallelCPUs)) tdm$parallelCPUs <- 1;
  if (is.null(tdm$path)) tdm$path <- paste(getwd(),"/",sep="");
  if (!is.null(tdm$runList)) {
    if (is.null(tdm$spotList)) tdm$spotList <- tdm$runList;
    if (is.null(tdm$filenameEnvT)) tdm$filenameEnvT=sub(".conf",".RData",tdm$runList[1],fixed=TRUE);
  }  
  
  # code which was previously in unbiasedRun. Now we put it here and call tdmDefaultsFill from unbiasedRun
  # (cleaner code, less places where tdm values are set)
  if (is.null(tdm$test2.string)) tdm$test2.string="default cutoff";
  if (is.null(tdm$tstCol)) tdm$tstCol="TST";
  if (!is.null(tdm$tstFrac)) {
    # this section only as long as tdm$tstFrac is still def'd in some script_* or other files
    cat("NOTE: Deprecated value tdm$tstFrac used. It will overwrite tdm$TST.testFrac.\n");
    tdm$TST.testFrac = tdm$tstFrac;
  }
  if (is.null(tdm$TST.testFrac)) tdm$TST.testFrac=0.2;
  if (is.null(tdm$nfold)) tdm$nfold=10;
  if (is.null(tdm$nrun)) tdm$nrun=5;
  if (is.null(tdm$U.saveModel)) tdm$U.saveModel=FALSE;

  tdm;
}
