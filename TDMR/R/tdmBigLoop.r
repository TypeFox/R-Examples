require(SPOT);

######################################################################################
# tdmBigLoop:
#
#' Tuning and unbiased evaluation in a big loop.
#'
#' For each \code{.conf} file in \code{tdm$runList} call all tuning algorithms (SPOT, CMA-ES or other) specified in \code{tdm$tuneMethod}
#' (via function \code{\link{tdmDispatchTuner}}). After each tuning process perform a run 
#' of \code{tdm$unbiasedFunc} (usually \code{\link{unbiasedRun}}). \cr
#' Each of these experiments is repeated tdm$nEperim times. Thus we have for each tripel \cr
#'    \tabular{ll}{
#'        \tab \code{   (confFile,nExp,theTuner)} \cr
#'    }
#' a tuning result. The ranges of the triple elements are: \cr
#'    \tabular{ll}{
#'        \tab \code{   confFile in tdm$runList} \cr
#'        \tab \code{   nExp in 1,...,tdm$nExperim} \cr
#'        \tab \code{   theTuner in tdm$tuneMethod} \cr
#'    }
#' 
#Details:
#' \code{tdm} refers to \code{envT$tdm}.
#' \cr\cr
#' Tuning is skipped if the \code{.conf} file does not appear in \code{tdm$spotList} or if \code{spotStep=="rep"}. In this
#' case it is assumed then that \code{envT$bstGrid} and \code{envT$resGrid} contain the appropriate data frames already.
#' See \code{\link{tdmEnvTAddBstRes}} on how to fill \code{envT$bstGrid} and \code{envT$resGrid}  from an \code{.RData} file.
#'
#' The available tuning algorithms (tuners) are 
#'      \itemize{
#'      \item{\code{\link{spotTuner}}:  Call \code{\link{spot}}.   }
#'      \item{\code{\link{lhdTuner}}:  Perform a parameter tuning using a Latin hypercube design (LHD) 
#'            for obtaining best design points. LHD is performed by configuring SPOT 
#'            in such a way that all the budget is used for the initial design (usually LHD). }
#'      \item{\code{\link{cma_jTuner}}:  Perform a parameter tuning by CMA-ES, using the *Java* 
#'            implementation by Niko Hansen through the interface package \code{\link[rCMA]{rCMA}}.    }
#'      \item{\code{\link{cmaesTuner}}:  Perform a parameter tuning by CMA-ES, using the *R*-implementation 
#'            (package \href{http://cran.r-project.org/web/packages/cmaes/}{\code{cmaes}} by Olaf Mersmann) 
#'            (deprecated, use \code{\link{cma_jTuner}} instead).  }
#'      \item{\code{\link{bfgsTuner}}:   Perform a parameter tuning by Broyden, Fletcher, Goldfarb and Shanno (BFGS) method.
#'            The L-BFGS-B version allowing box constraints is used.  }
#'      \item{\code{\link{powellTuner}}:  Perform a parameter tuning by Powell's UObyQA algorithm 
#'            (unconstrained optimization by quadratic approximation), see package \code{\link[powell]{powell}}).   } 
#'      }
#'
#'  @param envT      an environment containing on input at least the element \code{tdm} (a list with general settings for TDMR, 
#'                   see \code{\link{tdmDefaultsFill}}), which has at least the elements  
#'     \describe{
#'     \item{\code{tdm$runList}}{ vector of \code{.conf} filenames }
#'     \item{\code{tdm$spotList}}{ \code{[NULL]} vector of \code{.conf} filenames for which spot tuning is done. 
#'                             If \code{NULL}, then \code{spotList=runList} }
#'     }
#'  @param spotStep  \code{["auto"]} which step of SPOT to execute (either \code{"auto"} or \code{"rep"}).
#'  @param dataObj   \code{[NULL]} optional object of class \code{\link{TDMdata}} (the same for all runs in big loop). 
#'      If it is \code{NULL}, it will be constructed here with the help of \code{\link{tdmSplitTestData}}.
#'      Then it may be different for each run in the big loop.
#'
#'  @return environment \code{envT}, containing  the results
#'      \item{res}{ data frame with results from last tuning (one line for each call of \code{tdmStart*})} 
#'      \item{bst}{ data frame with the best-so-far results from last tuning (one line collected after each (SPO) step)}
#'      \item{resGrid}{  list with data frames \code{res} from all tuning runs. Use \cr
#'            \code{envT$getRes(envT,confFile,nExp,theTuner)}  \cr
#'        to retrieve a specific \code{res}. }
#'      \item{bstGrid}{  list with data frames \code{bst} from all tuning runs. Use \cr
#'            \code{envT$getBst(envT,confFile,nExp,theTuner)}  \cr
#'        to retrieve a specific \code{bst}. }
#'      \item{theFinals}{ data frame with one line for each triple \code{(confFile,nExp,tuner)}, each line contains summary
#'        information about the tuning run in the form: \cr
#'            \code{confFile tuner nExp [params] NRUN NEVAL RGain.bst RGain.* sdR.*} \cr
#'        where \code{[params]} is written depending on \code{tdm$withParams}. \cr
#'        \code{NRUN} is the number of unbiased evaluation runs. \cr
#'        \code{NEVAL} is the number of function evaluations (model builds) during tuning. \cr
#'        \code{RGain} denotes the relative gain on a certain data set: the actual gain achieved with the model 
#'        divided by the maximum gain possible for the current cost matrix and the current data set. This is for classification
#'        tasks, in the case of regression each \code{RGain.*} is replaced by \code{RMAE.*}, the relative mean absolute error. \cr
#'        Each 'sdR.' denotes the standard deviation of the preceeding RGain or RMAE. \cr
#'        RGain.bst is the best result during tuning obtained on the training-validation data. RGain.avg is the average result 
#'        during tuning. The following pairs {RGain.* sdR.*} are the results of one or several unbiased evaluations on the test data
#'        where '*' takes as many values as there are elements in \code{tdm$umode} (the possible values are explained in 
#'        \code{\link{unbiasedRun}}).  
#'        }
#'      \item{result}{ object of class \code{\link{TDMclassifier}} or \code{\link{TDMregressor}}. This is a list with results from \code{tdm$mainFunc} 
#'          as called in the last unbiased evaluation using the best parameters found during tuning. 
#'          Use \code{\link[=print.TDMclassifier]{print}(envT$result)} to get more info on such an object of class \code{\link{TDMclassifier}}.  }
#'      \item{tunerVal}{ an object with the return value from the last tuning process. For every tuner, this is the list 
#'          \code{spotConfig}, containing the SPOT settings plus the TDMR settings in elements \code{opts} and \code{tdm}. Every tuner 
#'          extends this list by \code{tunerVal$alg.currentResult} and \code{tunerVal$alg.currentBest}, see \code{\link{tdmDispatchTuner}}.
#'          In addition, each tuning method might add specific elements to the list, see the description of each tuner. }
#'   Environment \code{envT} contains further elements, but they are only relevant for the internal operation of 
#'   \code{tdmBigLoop} and its subfunctions.
#'
#' @note Side effects:
#'   Irrespective of the value of \code{tdm$fileMode}, 
#'     \itemize{
#'         \item a compressed version of \code{envT } is saved to file \code{tdm$filenameEnvT} (default: \code{<runList[1]>.RData}), 
#'               relative to  the directory of the \code{.conf} file. 
#'     }
#'   If \code{tdm$U.saveModel==TRUE}, then \code{envT$result$lastRes$lastModel} (the last trained model) will be saved to \code{tdm$filenameEnvT}. 
#'   The default is \code{tdm$U.saveModel==FALSE} (smaller \code{.RData} files).
#' 
#'   If \code{tdm$fileMode==TRUE}, more files are written relative to  the directory of the \code{.conf} file:
#'     \itemize{
#'         \item \code{envT$theFinals } is written to file \code{tdm$finalFile} 
#'         \item \code{envT$theFinals } is appended to \code{tdm$experFile}
#'      }                                                                                                           
#'   If \code{tdm$finalFile==NULL}, then it is set to \code{sub(".conf",".fin",runList[1])}.  \cr
#'   If \code{tdm$experFile==NULL}, then nothing is appended to any experiment file.
#'
#' Example usages of function \code{tdmBigLoop} are shown in \cr
#' \tabular{ll}{
#'    \tab \code{   demo(demo03sonar)} \cr 
#'    \tab \code{   demo(demo03sonar_B)} \cr
#'    \tab \code{   demo(demo04cpu)} \cr
#' }
#' where the corresponding R-sources are in directory \code{demo}.
#'
#' @examples
#' #*# This demo shows a complete tuned data mining process (level 3 of TDMR) where 
#' #*# the data mining task is the classification task SONAR (from UCI repository, 
#' #*# http://archive.ics.uci.edu/ml/datasets/Connectionist+Bench+%28Sonar,+Mines+vs.+Rocks%29).
#' #*# The data mining process is in main_sonar.r, which calls tdmClassifyLoop and tdmClassify
#' #*# with Random Forest as the prediction model. 
#' #*# The three parameter to be tuned are CUTOFF1, CLASSWT2 and XPERC, as specified 
#' #*# in file sonar_04.roi. The tuner used here is LHD.  
#' #*# Tuning runs are rather short, to make the example run quickly. 
#' #*# Do not expect good numeric results. 
#' #*# See demo/demo03sonar_B.r for a somewhat longer tuning run, with two tuners SPOT and LHD.
#' 
#' ## set working directory (dir with .apd, .conf and main_*.r file)
#' path <- paste(find.package("TDMR"), "demo02sonar",sep="/");
#' source(paste(path,"main_sonar.r",sep="/"));    
#' 
#' ## control settings for TDMR
#' tdm <- list( mainFunc="main_sonar"
#'            , runList = c("sonar_04.conf")
#'            , umode="CV"              # { "CV" | "RSUB" | "TST" | "SP_T" }
#'            , tuneMethod = c("lhd")
#'            , filenameEnvT="exBigLoop.RData"   # file to save environment envT (in dir 'path')
#'            , nrun=1, nfold=2         # repeats and CV-folds for the unbiased runs
#'            , nExperim=1
#'            , parallelCPUs=1
#'            , parallelFuncs=c("readCmdSonar")
#'            , optsVerbosity = 0       # the verbosity for the unbiased runs
#'            );
#' ## Each element of tdm$runList has the settings for one tuning process (e.g. 
#' ##    - auto.loop.steps = number of SPOT generations       
#' ##    - auto.loop.evals = budget of model building runs and 
#' ##    - io.roiFileName = "sonar_04.roi"
#' ## ). 
#' 
#' spotStep = "auto";   
#' source(paste(path,"start_bigLoop.r",sep="/"),chdir=TRUE);    # change dir to 'path' while sourcing
#'
#' @seealso   \code{\link{tdmDispatchTuner}}, \code{\link{unbiasedRun}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de}), Patrick Koch
#' @export
######################################################################################
tdmBigLoop <- function(envT,spotStep="auto",dataObj=NULL) {
                            # The environment envT is passed by reference into the inner functions
                            # which means that it can be used a) to transport information back from
                            # those inner functions and b) to transport information to and back even 
                            # for functions like tdmStartOther which are not allowed to have envT in 
                            # their argument list
  
  tdm <- envT$tdm;
  tdm <- tdmMapDesLoad(tdm);

  if (is.null(envT$getBst)) envT <- tdmEnvTAddGetters(envT);
  # just in case the user loaded envT via load("myFile.Rdata") and not tdmEnvTLoad("myFile.Rdata") 
      
  if (!is.null(tdm$mainFile)) source(tdm$mainFile);           # deprecated (sourcing main_TASK should be done in caller of tdmBigLoop)
  nTuner <- length(tdm$tuneMethod);
  nRunList <- length(tdm$runList);

	# depending on tdm$parallelCPUs, the following lines execute the code either in parallel 
	# or sequential: For each value of indVec the function bigLoopStep is 
	# called. The vectors tuneVec, expeVec, confVec contain for each element of indVec the
	# corresponding value of tuner, experiment number and .conf file, resp.
	# In case of parallel execution, parSapply will only return after the last parallel job has finished.
	#
	indVec <- 1:(tdm$nExperim*nTuner*nRunList);
	tuneVec <- rep(tdm$tuneMethod,tdm$nExperim*nRunList);
	expeVec <- rep(sort(rep(1:tdm$nExperim,nTuner)),nRunList);
	confVec <- sort(rep(tdm$runList,tdm$nExperim*nTuner));
	if (length(spotStep)>1) stop("Only a scalar string allowed for spotStep");
	if (length(tdm$umode)>1) stop("Only a scalar string allowed for tdm$umode");
	if (length(indVec)==1 & tdm$parallelCPUs>1) {
	  warning("There is only one job (length(indVec)==1) --> parallelization would not run, so we set tdm$ParallelCPUs=1");
	  tdm$parallelCPUs=1;
  }

  if (tdm$parallelCPUs>1) {
    cl = prepareParallelExec(tdm);
    sappResult <- parallel::parSapply(cl, indVec, bigLoopStep, tuneVec,expeVec,confVec,tdm$spotList,spotStep,dataObj,envT,tdm);
    parallel::stopCluster(cl);
    #--- old version was:
    #prepareParallelExec_SF(tdm);   # this is now in deprecated sources
    #sappResult <- sfSapply(indVec, fun=bigLoopStep, tuneVec,expeVec,confVec,tdm$spotList,spotStep,dataObj,envT,tdm);
  } else {  
		sappResult <- sapply(indVec, bigLoopStep, tuneVec,expeVec,confVec,tdm$spotList,spotStep,dataObj,envT,tdm);
  }
  # populate envT with the results returned in matrix sappResult:
  envT <- populateEnvT(sappResult,envT,tdm,spotStep);
  
  if (spotStep == "auto") {
    saveEnvT(envT,tdm$runList,tdm$filenameEnvT,saveModel=tdm$U.saveModel);
    saveSRFinfo(envT);
  } 

  if (nrow(envT$theFinals)>0) {
    print(envT$theFinals);
  } else {
    cat("Note: No rows in data frame envT$theFinals\n");
  }

  envT;
} # function tdmBigLoop

      ######################################################################################
  		#------------------------------------------------------------------------------------------------------
      #  bigLoopStep: helper function for tdmBigLoop, called via sapply or parSapply:
      #  (ind is an index where confFile varies slowest, nExp varies 2nd-slowest and theTuner varies fastest)
      bigLoopStep <- function(ind,tuneVec,expeVec,confVec,spotList,spotStep,dataObj,envT,tdm) {
        if (tdm$parallelCPUs>1) library(TDMR);
        theTuner = tuneVec[ind];
        nExp = expeVec[ind];
        confFile = confVec[ind];
        nConf = which(confFile==envT$runList);
        #i <- (nConf-1) %% length(spotStep) + 1; # i is an index which cyclically re-uses entries      
        #                                        # from vector spotStep, if it is shorter than runList
        #--- deprecated, now spotStep is only a string, "auto" or "rep" or "report"
        
        print(c(ind,confFile,nExp,theTuner));
        envT$spotConfig <- sC <- envT$sCList[[nConf]]; # spotGetOptions(srcPath=tdm$theSpotPath,confFile);
        envT$theTuner <- theTuner;
        envT$nExp <- nExp;
        envT$bst <- NULL;
        envT$res <- NULL; 
        if (spotStep=="rep" | spotStep=="report" | !(confFile %in% spotList)) {
          envT$bst = envT$getBst(envT,confFile,nExp,theTuner);
          envT$res = envT$getRes(envT,confFile,nExp,theTuner);
        }

        if (tdm$fileMode) {  
          if (is.null(tdm$finalFile)) tdm$finalFile = sub(".conf",".fin",confFile);
          tFinalFile <- ifelse(tdm$nExperim>1, sprintf("%s-e%02d%s",sub(".fin","",tdm$finalFile,fixed=TRUE),nExp,".fin"), tdm$finalFile);
          # i.e. if tdm$finalFile="cpu.fin", then tFinalFile="cpu-e02.fin" for nExp=2
          if (file.exists(tFinalFile) & theTuner==tuneVec[1] & confFile==envT$runList[1]) file.remove(tFinalFile);
        } 
        # NEW 09/2012: always operate SPOT with spot.fileMode=FALSE (take everything from  envT$spotConfig)
        envT$spotConfig$spot.fileMode=FALSE;
        
        #
        # this is the preferred place to read the data and split them into test data and train/vali data
        #
        if (is.null(dataObj)) dataObj <- tdmSplitTestData(sC$opts,tdm,nExp);
        if (!is.null(dataObj)) envT$spotConfig$opts$TST.COL = dataObj$TST.COL;    # this column has to be subtracted in main_* from the input variables

        ptm <- proc.time();
        if (confFile %in% spotList) {
            # this is for the case spotStep=="rep":
            envT$spotConfig$alg.currentResult <- envT$res;
            envT$spotConfig$alg.currentBest <- envT$bst;	
            cat(sprintf("*** Starting TUNER %s, spotStep=%s, on task %s ***\n",theTuner,spotStep,confFile));
            tdmDispatchTuner(theTuner,confFile,spotStep,tdm,envT,dataObj);
            # If spotStep=="auto" then tdmDispatchTuner runs the tuning process, puts its results in envT$bst, envT$res 
            # and returns the (extended) spotConfig in envT$tunerVal.
            # If spotStep=="rep" then tdmDispatchTuner expects that envT$bst, envT$res contain already the right data frames 
            # and it runs just the reporting step of SPOT.
        } 
        else {    # i.e. !(confFile %in% spotList)
            # we assume that envT$bst and envT$res contain already the right data frames (filled in from a previous setting)
        }    
        
        #
        # now we have envT$bst and envT$res filled in any case (either from tdmDispatchTuner or from a previous setting)
        #
        if (is.null(envT$bst)) stop("No contents for data frame envT$bst");
        if (is.null(envT$res)) stop("No contents for data frame envT$res");
        
        if (is.null(tdm$timeMode)) stop("tdm$timeMode is not set (NULL). Consider 'tdm <- tdmDefaultsFill(tdm)' to set all defaults");
        time.txt = c("Proc", "System", "Elapsed");
        time.TRN=(proc.time()-ptm)[tdm$timeMode]; opts2=list(); opts2$VERBOSE=1;  
        cat1(opts2,paste(time.txt[tdm$timeMode], "time for tuning with tdmDispatchTuner:",time.TRN,"sec\n"));     

        ptm <- proc.time();
        finals <- NULL;
        if (tdm$nrun>0) {
          for (umode in tdm$umode) {
            cat("*** starting",tdm$unbiasedFunc,"for",confFile,"with umode=",umode,"***\n");
            cmd=paste("finals <-",tdm$unbiasedFunc,"(confFile,envT,dataObj,umode=umode,finals=finals,withParams=envT$wP,tdm=tdm)",sep="");
            eval(parse(text=cmd));
          }
          time.TST=(proc.time()-ptm)[tdm$timeMode];   
          cat1(opts2,paste(time.txt[tdm$timeMode], "time for",tdm$unbiasedFunc,":",time.TST,"sec\n"));   
          finals <- cbind(finals,Time.TST=time.TST,Time.TRN=time.TRN);  

          # transport back opts$srf (data frame with SRF info in case opts$SRF.calc==TRUE)
          envT$sCList[[nConf]]$opts$srf = envT$result$lastRes$opts$srf

          if (tdm$fileMode) {  
            #removeTmpfiles2(confFile);
            if (!file.exists(dirname(tFinalFile))) {
              success = dir.create(dirname(tFinalFile));     
              if (!success) stop(sprintf("Could not create dirname(tFinalFile)=%s",dirname(tFinalFile)));
            }
        		colNames = ifelse(file.exists(tFinalFile),FALSE,TRUE);
        		write.table(finals
        				, file = tFinalFile
        				, col.names= colNames
        				, row.names= FALSE
        				, append = !colNames
        				, sep = " ",
        				, quote = FALSE
        				, eol = "\n"
        		);
        		if (!is.null(tdm$experFile)) {
              if (!file.exists(dirname(tdm$experFile))) {
                success = dir.create(dirname(tdm$experFile));     
                if (!success) stop(sprintf("Could not create dirname(tdm$experFile)=%s",dirname(tdm$experFile)));
              }
          		colNames = ifelse(file.exists(tdm$experFile),FALSE,TRUE);
          		write.table(finals             # multiple execution of same experiment (nExperim>1 in script_all.R)
          				, file = tdm$experFile
          				, col.names= colNames
          				, row.names= FALSE
          				, append = !colNames
          				, sep = " ",
          				, quote = FALSE
          				, eol = "\n"
          		);
        		}
          } # if(tdm$fileMode)
          flush.console(); 
        } # if(tdm$nrun>0)
         
        # we return a *list* with the most important elements from environment envT and not envT itself, because 
        # parSapply (parallel execution) can not deal with environments as return values.
        list( tunerVal=envT$tunerVal    # last tuning result: spotConfig with extensions       
             ,bst=envT$bst              # last tuning result       
             ,res=envT$res              # last tuning result  
             ,result=envT$result        # last tuning result
             ,roi=envT$tunerVal$alg.roi      
             ,theFinals=finals   ###tail(envT$theFinals,1)
             );      
  		} # end of function bigLoopStep
  		#------------------------------------------------------------------------------------------------------

######################################################################################
# helper fct for tdmBigLoop: 
#     Populate the global envT after parallel execution  with 'parallel' (parSapply with results in sappResult).
#     For simplicity we use it also after sequential execution (sapply with results in sappResult). 
populateEnvT <- function(sappResult,envT,tdm,spotStep) {
    if (spotStep=="auto") {
      nGrid = length(envT$bstGrid);
      for (ind in 1:(tdm$nExperim*length(tdm$tuneMethod)*length(envT$runList))) {
          nGrid = nGrid+1;
          envT$bstGrid[[nGrid]] <- as.data.frame(sappResult["bst",ind][[1]]);
          envT$resGrid[[nGrid]] <- as.data.frame(sappResult["res",ind][[1]]);   
          envT$roiGrid[[nGrid]] <- as.data.frame(sappResult["roi",ind][[1]]);   
          envT$theFinals <- rbind(envT$theFinals,as.data.frame(sappResult["theFinals",ind][[1]]));
          # "[[1]]" is necessary to avoid prefix "theFinals." in the header names 
      }
      if (!is.null(envT$theFinals)) if(nrow(envT$theFinals)>0)
          rownames(envT$theFinals) <- (1:nrow(envT$theFinals));
      envT$bst <- as.data.frame(sappResult["bst",ncol(sappResult)][[1]]);
      envT$res <- as.data.frame(sappResult["res",ncol(sappResult)][[1]]);
      envT$tunerVal <- sappResult["tunerVal",ncol(sappResult)][[1]];    # last tuning result
      envT$result <- sappResult["result",ncol(sappResult)][[1]];        # last tuning result
    }
    else {     # i.e. if spotStep == "rep" or == "report"
      nGrid = nrow(envT$theFinals);
      if (ncol(sappResult) > nGrid) stop("There are more columns in sappResult than rows in envT$theFinals");
      for (ind in 1:ncol(sappResult)) {
        theFinals <- sappResult["theFinals",ind][[1]]
        names_finals <-  paste(c(Opts(envT$result)$rgain.string,"sdR"),tdm$umode,sep=".");
        envT$theFinals[ind,names_finals] <- theFinals[names_finals];
        envT$theFinals[ind,"Time.TST"] <- theFinals["Time.TST"];
      }
    }
    envT;
}

######################################################################################
# helper fct for tdmBigLoop: 
#     Save a small version of environment envT (passed in as thisEnvT) on filenameEnvT
#     (if NULL, use <runList[1]>.RData).
#     If savePredictions is TRUE, the elements envT$result$predictions, envT$result$predProbList, envT$result$lastRes$predProb are
#     saved. The default is savePredictions==FALSE.
#     If saveModel is TRUE, the element envT$result$lastRes$lastModel is saved. The default is savePredictions==FALSE.
saveEnvT <- function(thisEnvT,runList,filenameEnvT=NULL,savePredictions=FALSE,saveModel=FALSE) {
      envT = list() # new.env(); #            # when saving thisEnvT, we copy the relevant elements to a *list* envT
                                              # and skip the function elements (getBst,getInd,getRes) in envT, because 
                                              # they would also save *their* environment which tends to make the .RData file rather big
      class(envT) = c("TDMenvir")
      for (ele in setdiff(ls(thisEnvT),c("getBst","getInd","getRes"))) 
        eval(parse(text=paste("envT$",ele," = thisEnvT$",ele,sep="")));
                                              # for the save on .RData file we delete also some potentially voluminous elements 
      envT$result$lastRes$d_train=NULL;       # from the local envT, we want only res and bst objects + describing vars ...
      envT$result$lastRes$d_test=NULL;
      envT$result$lastRes$d_dis=NULL;
      envT$result$lastRes$lastProbs=NULL;
      if (!saveModel) envT$result$lastRes$lastModel=NULL;
      envT$result$dset=NULL;
      envT$result$TST=NULL;
      if (!savePredictions) {
        envT$result$predictions=NULL;     
        envT$result$predProbList=NULL;  
        envT$result$lastRes$predProb=NULL;   
      }
      envT$tunerVal$seq.modelFit=NULL;        # this can be quite big in the case of spotPredictRandomForest
      envT$spotConfig=NULL;     # versions of spotConfig are contained in envT$sCList[[i]] and envT$tunerVal
      envT$theTuner=NULL;
      envT$nExp=NULL;           # see envT$tdm$nExperim for number of experiments
      if (is.null(filenameEnvT)) filenameEnvT=sub(".conf",".RData",runList[1],fixed=TRUE);

      save(envT,file=filenameEnvT);              
                                              # ... but we leave thisEnvT (which is envT in tdmBigLoop) untouched
                                              #     and thus return the full envT environment to the caller of tdmBigLoop
}

######################################################################################
# helper fct for tdmBigLoop:
#     Save the info for SRF (sorted Random Forest) importance on file(s) in case opts$SRF.calc==TRUE
#     (Nothing is done in case opts$SRF.calc==FALSE.)
saveSRFinfo <- function(envT) {
  # function saveSRF saves the field opts$srf (a list with as many data frames as there are response variables)
  # in case opts$SRF.calc==TRUE to RData file opts$srfFile (filename set by function addSRF in tdmEnvTMakeNew.r)
  saveSRF <- function(opts) {
      if (opts$SRF.calc==TRUE) {
        cat("Saving sorted RF importance info to file", opts$srfFile, ".\n")
        srf = opts$srf;
        save(file=opts$srfFile,srf);
      }
  }
  k=0;
  for (confFile in envT$runList) {
    k=k+1;
 	  saveSRF(envT$sCList[[k]]$opts);
  }
}

######################################################################################
# helper fct for tdmBigLoop 
removeTmpfiles2 <- function(confFile) {
      #for (suf in c("apd","aroi","bst","conf","des","res","roi"))  {    # old version (removeTmpfiles)
      for (suf in c("aroi","bst","des","res"))  {
        tFile <- paste(strsplit(confFile,"conf"),suf,sep="");
        if (file.exists(tFile)) file.remove(tFile);
      }
}
   
######################################################################################
# helper function for tdmBigLoop, only needed in case tdmParallelCPUs>1 
# 
# set up everything for parallel execution
#
prepareParallelExec <- function(tdm) 
{
    # require(parallel);     # now via direct call 'parallel::' 

    cl <- parallel::makeCluster(tdm$parallelCPUs)
    parallel::clusterExport(cl,c("tdm"))
    
    if (!exists(".Random.seed")) set.seed(42);
    parallel::clusterExport(cl,".Random.seed");
    parallel::clusterExport(cl,c(tdm$mainFunc,tdm$parallelFuncs));
    if (is.null(tdm$tdmPath)) {
        cat("Export installed library TDMR to parallel cluster\n");
        #library(TDMR);     # this is now done in function bigLoopStep
    } else {
        stop("For tdm$parallelCPUs>1 it is required to use the *library* version of TDMR. Consider to set tdm$tdmPath=NULL.")
    }
    cl;
} # prepareParallelExec()
  

