require(SPOT);
#-------- MIES part temporarily moved to ../tdmDispatch_miesPart.r

###########################################################################################
#tdmDispatchTuner:
#
#'     Helper function for \code{\link{tdmBigLoop}}.
#'
#'     tdmDispatchTuner selects and starts the tuner specified by tuneMethod. \cr
#'     See the 'Details' section of \code{\link{tdmBigLoop}} for a list of available tuners.
#'
#' @param tuneMethod the tuning algorithm given as a string. Possible values are \{ "spot" | "lhd" | "cmaes" | "cma_j" | "bfgs" | "powell" \}.
#' @param confFile the configuration file.
#' @param spotStep which step to execute for \link{SPOT} . Values "rep" ("report") and "auto" are supported by TDMR.
#' @param tdm the TDMR object
#' @param envT the environment variable
#' @param dataObj the \code{\link{TDMdata}} object containing the data set (train/vali part and test part)
#'
#' @return The result \code{tunerVal} of the tuning algorithm, this is the list envT$spotConfig, extended by
#'      \item{\code{alg.currentResult}}{ the RES data frame  }
#'      \item{\code{alg.currentBest}}{ the BST data frame  }
#'      \item{...}{[optional] further results from the tuning algorithm \code{tuneMethod}  }
#'
#' @note Side effects: Environment \code{envT} is extended by
#'     \itemize{
#'         \item \code{envT$tunerVal = tunerVal} 
#'         \item \code{envT$res  } the RES data frame
#'         \item \code{envT$bst  } the BST data frame
#'      }
#' @seealso   \code{\link{tdmBigLoop}}
#' @author Wolfgang Konen, FHK, 2010 - 2013
#' @export
#' @keywords internal
###########################################################################################
tdmDispatchTuner <- function(tuneMethod,confFile,spotStep,tdm,envT,dataObj)
{
    theSpotPath=tdm$theSpotPath
    if (spotStep=="auto") {
      #if (is.null(tdm$mainFile)) stop("Element tdm$mainFile is missing, but this is required for spotStep='auto'");
      if (is.null(tdm$mainFunc)) stop("Element tdm$mainFunc is missing, but this is required for spotStep='auto'");
    }   
    if (!is.null(envT$spotConfig$alg.currentBest)) {
      # Fix for upgading older BST data frames to the new SPOT version (V1.0.2662): 
      # A BST data frame now has to have a column STEP, otherwise a crash in spotWriteBest (spotHelpFunctions.R) will occur.  
      if (!any(names(envT$spotConfig$alg.currentBest)=="STEP")) {
        cat("NOTE: Since envT$spotConfig$alg.currentBest has no column STEP, we add for SPOT a dummy column STEP containing the numbers 1 :",nrow(envT$spotConfig$alg.currentBest),"\n");
        envT$spotConfig$alg.currentBest <- cbind(envT$spotConfig$alg.currentBest,STEP=(1:nrow(envT$spotConfig$alg.currentBest)));
      }
    }
    if (tdm$fileMode==FALSE) confFile="NULL";      # everything is passed to spot() via envT$spotConfig
    
    tunerVal = switch(spotStep
      ,"rep"=,"report"= spot(confFile,spotStep,theSpotPath,envT$spotConfig)
      ,"auto"= switch(tuneMethod
                ,"spot" = spotTuner(confFile,spotStep,theSpotPath,tdm,envT,dataObj)
                ,"lhd" = lhdTuner(confFile,spotStep,theSpotPath,tdm,envT,dataObj) 
                ,"cmaes" = cmaesTuner(confFile,tdm,envT,dataObj)
                #,"mies" = miesTuner(confFile,tdm,envT,dataObj)  
                ,"cma_j" = cma_jTuner(confFile,tdm,envT,dataObj) 
                ,"bfgs" = bfgsTuner(confFile,tdm,envT,dataObj) 
                ,"powell" = powellTuner(confFile,tdm,envT,dataObj) 
                ,"INVALID1"
                )
      ,"INVALID2"
      );

    if (tunerVal[1]=="INVALID1") {
      stop(sprintf("*** Invalid tuneMethod=%s ***\n",tuneMethod));
    }
    if (tunerVal[1]=="INVALID2") {
      stop(sprintf("*** Invalid spotStep=%s ***\n",spotStep));
    }
    tunerVal$dataObj = NULL;     # delete this potential voluminous element

    if (!(tuneMethod %in% c("spot","lhd"))) {
      # Bug fix 05/12:        
      # When tuning via tdmStartOther terminates, there might have been last calls to tdmStartOther which have not
      # yet been processed for the BST data frame (or this data frame may not have been constructed yet).
      # So we compute again, based on tunerVal$alg.currentResult, the so far best solution (merge repeats), 
      # and write a line to the BST data frame tunerVal$alg.currentBest.
    	mergedData <- spotPrepareData(tunerVal)
      tunerVal <- spotWriteBest(mergedData, tunerVal);	# appends the best solution to tunerVal$alg.currentBest
  		envT$bst <- tunerVal$alg.currentBest;  
    }

    envT$tunerVal = tunerVal;
    tunerVal;
}

######################################################################################
#' Perform a \link{SPOT} tuning. 
#'
#' \code{spotTuner} calls \code{\link{spot}} and write the necessary information on environment \code{envT}.
#'
#' @param confFile task configuration for tuning algorithm
#' @param spotStep needed for spot
#' @param theSpotPath needed for spot
#' @param tdm the TDMR object
#' @param envT the environment variable
#' @param dataObj the \code{\link{TDMdata}} object containing the data set (train/vali part and test part)
#'
#' @return the result of SPOT tuning, i.e. the list \code{spotConfig}
#'
# @rdname  tdmTuners
# @aliases lhdTuner spotTuner
#'
#' @export
#' @keywords internal
######################################################################################
spotTuner <- function(confFile,spotStep,theSpotPath,tdm,envT,dataObj)
{
    sC <- envT$spotConfig;
    #sC$alg.func <- "tdmStartSpot";
    #sC$spot.fileMode = TRUE; 
    sC$tdm <- tdm;              # needed for tdm$mainFile, tdm$mainFunc, tdm$fileMode
    sC$alg.currentResult <- NULL;
    sC$dataObj <- dataObj;

    sC <- spot(confFile,spotStep,theSpotPath,sC);
    # spot calls tdmStartSpot. 
    # tdmStartSpot reads opts from sC$opts, dataObj from sC$dataObj (it takes only the train-validation part), 
    # and writes sC$alg.currentResult

    envT$res <- sC$alg.currentResult;
    envT$bst <- sC$alg.currentBest;	
    sC;
}

######################################################################################
#' Perform LHD tuning. 
#'
#' Perform a parameter tuning using a Latin hypercube design (LHD) 
#' for obtaining best design points. LHD is performed by configuring SPOT 
#' in such a way that all the budget is used for the initial design (usually LHD).
#' 
#' @param confFile task configuration for tuning algorithm
#' @param spotStep needed for spot
#' @param theSpotPath needed for spot
#' @param tdm the TDMR object
#' @param envT the environment variable
#' @param dataObj the \code{\link{TDMdata}} object containing the data set (train/vali part and test part)
#'
#' @return the result of LHD sampling, i.e. the list \code{spotConfig}
#'
#' @export
#' @keywords internal
######################################################################################
lhdTuner <- function(confFile,spotStep,theSpotPath,tdm,envT,dataObj)
{
    #if (tdm$fileMode==FALSE) {
    #  envT$spotConfig$spot.fileMode=FALSE;
    #  confFile="NULL";
    #}  

    sC <- envT$spotConfig;
    sC$tdm <- tdm;              # needed for tdm$mainFile, tdm$mainFunc, tdm$fileMode
    sC$alg.currentResult <- NULL;

    fncall1 = sC$auto.loop.nevals;
    fncall2 = sC$auto.loop.steps* sC$seq.design.new.size* sC$seq.design.maxRepeats +
              sC$init.design.size* sC$init.design.repeats  +
             (sC$seq.design.maxRepeats- sC$init.design.repeats);
    new.nevals = min(fncall1,fncall2); 
    #sC$alg.func = "tdmStartSpot";
    sC$auto.loop.steps = 0;
    sC$init.design.size = round(new.nevals/sC$seq.design.maxRepeats);
    sC$init.design.repeats = sC$seq.design.maxRepeats;
    sC$dataObj <- dataObj;

    sC <- spot(confFile,spotStep,theSpotPath,sC);
    # spot calls tdmStartSpot. 
    # tdmStartSpot reads opts from sC$opts, dataObj from sC$dataObj, and writes sC$alg.currentResult

    envT$res <- sC$alg.currentResult;
    envT$bst <- sC$alg.currentBest;	
    sC;
}


######################################################################################
#' Perform CMA-ES tuning (R version). 
#'
#' Perform a parameter tuning by CMA-ES, using the *R*-implementation 
#' (package \code{cmaes} by Olaf Mersmann).
#' 
#' @param confFile task configuration for tuning algorithm
#' @param tdm the TDMR object
#' @param envT the environment variable
#' @param dataObj the \code{\link{TDMdata}} object containing the data set (train/vali part and test part)
#'
#' @return the result of CMA-ES tuning, i.e. the list \code{envT$spotConfig}, extended by
#'    \item{\code{cma}}{ the return value from \code{\link[cmaes]{cma_es}}   }
#'    \item{\code{cma$count}}{ the number of calls to \code{tdmStartOther}  }
#'
#' @note 
#'    This tuner \code{cmaesTuner} is deprecated. Use instead \code{\link{cma_jTuner}}.
#' @export
#' @keywords internal
######################################################################################
cmaesTuner <- function(confFile,tdm,envT,dataObj)
{
    # require(cmaes)      # now via direct call 'cmaes::' 
    envT$spotConfig$alg.currentResult <- NULL;    
    sC <- envT$spotConfig;
    roiLower <- sC$alg.roi[,1];
    roiUpper <- sC$alg.roi[,2];
    set.seed(sC$spot.seed);
    param <- runif(length(roiLower),roiLower,roiUpper);#/2;    # start configuration
    if (tdm$fileMode) {
      tdm$resFile <-  sC$io.resFileName;   
      tdm$bstFile <-  sC$io.bstFileName;   
      if (file.exists(sC$io.resFileName)) file.remove(sC$io.resFileName);
      if (file.exists(sC$io.bstFileName)) file.remove(sC$io.bstFileName);
    }

    fitFunc <- function(x, opts=NULL) {
      ## tdmStartOther returns the fitness (mean(yres)) and saves other
      ## diagnostic results in envT$spotConfig, envT$res, envT$bst
      tdmStartOther(x,tdm,envT,dataObj,opts);
    }
    #--- old, deprecated: 
    #if(is.null(tdm$startOtherFunc)) tdmStartOther = makeTdmStartOther(tdm,envT,dataObj) else
    #  tdmStartOther = tdm$startOtherFunc;

    fncall1 = sC$auto.loop.nevals;
    fncall2 = sC$auto.loop.steps* sC$seq.design.new.size* sC$seq.design.maxRepeats 
            + sC$init.design.size* sC$init.design.repeats;
            + (sC$seq.design.maxRepeats- sC$init.design.repeats);
    control = list(mu=4,lambda=4);
    control$maxit=min(fncall1,fncall2)/sC$seq.design.maxRepeats/control$lambda;
    # Be aware, that control$lambda (# of offsprings) controls the # of fct calls, 
    # NOT control$mu. The number of calls to tdmStartOther will be N=ceiling(control$maxit)*control$lambda, 
    # the number of calls to the DM training fct main_TASK will be N*sC$seq.design.maxRepeats.
    control$maxit=max(round(control$maxit),1);     # BUG FIX 05/12: if control$maxit were NOT an integer, we would get 
            # from cma_es: "Error in eigen.log[iter, ] <- rev(sort(e$values)): subscript out of bounds"
    control$diag.pop=FALSE;
    control$diag.sigma=TRUE;
    control$diag.eigen=TRUE;

    cma <- cmaes::cma_es(param,fitFunc,sC$opts,lower=roiLower,upper=roiUpper,control=control);
    # cma_es calls tdmStartOther and tdmStartOther appends to data frames envT$res and envT$bst
    # (and to envT$spotConfig$alg.currentResult and ...$alg.currentBest as well).
    # If tdm$fileMode==TRUE, then tdmStartOther will write envT$res to envT$spotConfig$io.resFileName
    # and envT$bst to envT$spotConfig$io.bstFileName 
    
    tunerVal = envT$spotConfig;
    tunerVal$cma = cma;
    tunerVal$cma$count=nrow(envT$res);
    
    cat(sprintf("Function calls to tdmStartOther: %d\n",tunerVal$cma$count[1]));

    tunerVal;
}


######################################################################################
# cma_jTuner
#
#' Perform CMA-ES tuning (Java version). 
#' 
#' Perform a parameter tuning by CMA-ES, using the *Java* 
#' implementation by Niko Hansen through the interface package \code{\link[rCMA]{rCMA}}.
#' 
#' @param confFile task configuration for tuning algorithm
#' @param tdm the list with TDMR settings. \cr
#'    If \code{tdm$CMA.propertyFile} is not NULL, then \code{\link[rCMA]{rCMA}} will read this CMA property file. If NULL, 
#'    the default file CMAEvolutionStrategy.properties from find.package("rCMA") is taken. \cr
#'    If \code{tdm$CMA.populationSize} is not NULL, the CMA population size will be set accordingly. If NULL, 
#'    the default population size \code{4 + 3*log(N)} is taken with \code{N} = number of parameters in .roi file.
#' @param envT the environment variable
#' @param dataObj the \code{\link{TDMdata}} object containing the data set (train/vali part and test part)
#'
#' @return the result of CMA-ES tuning, i.e. the list \code{envT$spotConfig}, extended by list \code{cma} with elements
#'    \item{\code{cma$sres}}{ a string vector with all console output from the optimization run \code{cmaOptimDP} of package \code{\link[rCMA]{rCMA}}  }
#'    \item{\code{cma$count}}{ the number of calls to \code{tdmStartOther}  }
#' @export
#' @keywords internal
######################################################################################
cma_jTuner <- function(confFile,tdm,envT,dataObj)
{
    #require(rCMA);         # now via direct call 'rCMA::' 
    tunerVal <- cma_jInternRCma(tdm,envT,dataObj);
    #tunerVal <- cma_jInternRJava(tdm,envT,dataObj);    # -- deprecated, now in R-DM-Template-deprecated\
    #tunerVal <- cma_jIntern__Old(tdm,envT,dataObj);    # -- deprecated, now in R-DM-Template-deprecated\
}
# 
# Internal function for the usage of CMA-ES (Java) from R (based on package rCMA)
#
cma_jInternRCma <- function(tdm,envT,dataObj) {
    #require(rCMA);         # now via direct call 'rCMA::'
    envT$spotConfig$alg.currentResult <- NULL;
    sC <- envT$spotConfig;
    roiLower <- sC$alg.roi[,1];
    roiUpper <- sC$alg.roi[,2];
    if (!(length(roiLower)==length(roiUpper))) stop("Length of roiLower and roiUpper differ");
    set.seed(sC$spot.seed);
    initialX <- runif(length(roiLower),roiLower,roiUpper);#/2;    # initial point (mean for population)
    stdDevs <- (roiUpper-roiLower)/10;                            # initial standard deviations
    
    ## we assume that sC$auto.loop.nevals is the limiting factor in SPOT tuning
    stopMaxFunEvals = sC$auto.loop.nevals/sC$seq.design.maxRepeats;           

    fitFunc <- function(x) {
        ## tdmStartCma_j (see tdmStartOther.r) returns the fitness (mean(yres)) and 
        ## saves other diagnostic results in envT$spotConfig, envT$res, envT$bst
        tdmStartCma_j(x,tdm,envT,dataObj);
    }
    isFeasible <- function(x) {
        if (any(x<roiLower) | any(x>roiUpper)) return(FALSE);
        return(TRUE);
    }
    
    cma <- rCMA::cmaNew(propFile=tdm$CMA.propertyFile);
    # if tdm$CMA.propertyFile is NULL, the default file CMAEvolutionStrategy.properties from find.package("rCMA") will be loaded
    
    if (!is.null(tdm$CMA.populationSize)) {
      rCMA::cmaSetPopulationSize(cma,tdm$CMA.populationSize);
    }
    rCMA::cmaSetStopMaxFunEvals(cma,stopMaxFunEvals);
    rCMA::cmaSetStopFitness(cma,-1e300);  # never stop due to fitness function value
    rCMA::cmaInit(cma,seed=sC$spot.seed,dimension=length(initialX),initialX=initialX, initialStandardDeviations=stdDevs);
    res1 = rCMA::cmaOptimDP(cma,fitFunc,isFeasible,verbose=2,iterPrint=1);
    bestSolution=rCMA::cmaEvalMeanX(cma,fitFunc,isFeasible);

    tunerVal = envT$spotConfig;
		tunerVal$alg.currentResult <- envT$res;
		tunerVal$alg.currentBest <- envT$bst;
		
    tunerVal$cma = list()
    tunerVal$cma$sres <- res1$sMsg; 
    tunerVal$cma$count=nrow(envT$res);
    cat(sprintf("Function calls to tdmStartOther: %d\n",tunerVal$cma$count[1]));

    tunerVal;
}


######################################################################################
#' Perform Powell's tuning. 
#'
#' Perform a parameter tuning by Powell's UObyQA algorithm 
#' (unconstrained optimization by quadratic approximation), see R-package \code{powell}
#' 
#' @param confFile task configuration for tuning algorithm
#' @param tdm the TDMR object
#' @param envT the environment variable
#' @param dataObj the \code{\link{TDMdata}} object containing the data set (train/vali part and test part)
#'
#' @return the result of Powell tuning, i.e. the list \code{envT$spotConfig}, extended by
#'    \item{\code{powell}}{ the return value from \code{\link[powell]{powell}}   }
#'    \item{\code{powell$count}}{ the number of calls to \code{tdmStartOther}  }
#' @export
#' @keywords internal
######################################################################################
powellTuner <- function(confFile,tdm,envT,dataObj){
    # require(powell)       # now via direct call 'powell::'
  
    #function to fix constraint violating steps:
    tdm$constraintFnc <- function(x,tdm) {
      for(i in 1:length(x)){
        if(x[i] < envT$spotConfig$alg.roi[i,1]) {x[i] <- envT$spotConfig$alg.roi[i,1]}
        if(x[i] > envT$spotConfig$alg.roi[i,2]) {x[i] <- envT$spotConfig$alg.roi[i,2]}
      }
      x;
    }
    
    #if (tdm$fileMode==FALSE) envT$spotConfig$spot.fileMode=FALSE;
    envT$spotConfig$alg.currentResult <- NULL;
 
    sC <- envT$spotConfig;
    roiLower <- sC$alg.roi[,1];
    roiUpper <- sC$alg.roi[,2];
  
    if (tdm$fileMode) {
      tdm$resFile <-  sC$io.resFileName;   
      tdm$bstFile <-  sC$io.bstFileName;   
      if (file.exists(sC$io.resFileName)) file.remove(sC$io.resFileName);
      if (file.exists(sC$io.bstFileName)) file.remove(sC$io.bstFileName);
    }
  
    fitFunc <- function(x,opts=NULL) {
      ## tdmStartOther returns the fitness (mean(yres)) and saves other
      ## diagnostic results in envT$spotConfig, envT$res, envT$bst
      tdmStartOther(x,tdm,envT,dataObj,opts);
    }
    #--- old, deprecated: 
    #if(is.null(tdm$startOtherFunc)) tdmStartOther = makeTdmStartOther(tdm,envT,dataObj) else
    #  tdmStartOther = tdm$startOtherFunc;


    param <- (roiLower+roiUpper)/2;    # start configuration
  
    maxEvaluations = sC$auto.loop.nevals/sC$seq.design.maxRepeats;

    resPowell <- powell::powell(param, fn=fitFunc, control=list(maxit=maxEvaluations), check.hessian=FALSE, opts=sC$opts);
    # powell calls tdmStartOther and tdmStartOther writes envT$res and envT$bst
    # (and to envT$spotConfig$alg.currentResult and ...$alg.currentBest as well).

    tunerVal = envT$spotConfig;
		tunerVal$alg.currentResult <- envT$res;       
		tunerVal$alg.currentBest <- envT$bst;         
    tunerVal$powell = resPowell;
    tunerVal$powell$count=nrow(envT$res);
    tunerVal;  
}


######################################################################################
#' Perform BFGS tuning. 
#'
#' Perform a parameter tuning by Broyden, Fletcher, Goldfarb and Shanno (BFGS) method.
#' The L-BFGS-B version allowing box constraints is used.
#' 
#' @param confFile task configuration for tuning algorithm
#' @param tdm the TDMR object
#' @param envT the environment variable
#' @param dataObj the \code{\link{TDMdata}} object containing the data set (train/vali part and test part)
#'
#' @return the result of BFGS tuning, i.e. the list \code{envT$spotConfig}, extended by
#'    \item{\code{bfgs}}{ the return value from \code{optim(...,method="L-BFGS-B")}   }
#'    \item{\code{bfgs$count}}{ the number of calls to \code{tdmStartOther}  }
#' @author Wolfgang Konen \email{wolfgang.konen@@fh-koeln.de}, Patrick Koch 
#' @export
#' @keywords internal
######################################################################################
bfgsTuner <- function(confFile,tdm,envT,dataObj){
    #require("optimx")
    envT$spotConfig$alg.currentResult <- NULL;
    #if (tdm$fileMode==FALSE) envT$spotConfig$spot.fileMode=FALSE;
    sC <- envT$spotConfig;
    roiLower <- sC$alg.roi[,1];
    roiUpper <- sC$alg.roi[,2];
    param <- runif(length(roiLower), min=roiLower, max=roiUpper) # create start vector generating uniformly randomized vector
     
    if (tdm$fileMode) {
      tdm$resFile <-  sC$io.resFileName;   
      tdm$bstFile <-  sC$io.bstFileName;   
      if (file.exists(sC$io.resFileName)) file.remove(sC$io.resFileName);
      if (file.exists(sC$io.bstFileName)) file.remove(sC$io.bstFileName);
    }
  
    fitFunc <- function(x, opts=NULL) {
      ## tdmStartOther returns the fitness (mean(yres)) and saves other
      ## diagnostic results in envT$spotConfig, envT$res, envT$bst
      tdmStartOther(x,tdm,envT,dataObj,opts);
    }
    #--- old, deprecated: 
    #if(is.null(tdm$startOtherFunc)) tdmStartOther = makeTdmStartOther(tdm,envT,dataObj) else
    #  tdmStartOther = tdm$startOtherFunc;

    tdm$roi <- sC$alg.roi;
  
    maxEvaluations = sC$auto.loop.nevals/sC$seq.design.maxRepeats;
    maxEvaluations = round(maxEvaluations);

    bfgs <- optim(par=param, fn=fitFunc, gr=NULL, method="L-BFGS-B", lower=roiLower, upper=roiUpper, control=list(maxit=maxEvaluations), opts=sC$opts)
    # optim calls tdmStartOther and tdmStartOther writes envT$res and envT$bst
    # (and to envT$spotConfig$alg.currentResult and ...$alg.currentBest as well).
  
    cat(">>>>>> RESULT:\n")
    print(bfgs)
  
    tunerVal = envT$spotConfig;
		tunerVal$alg.currentResult <- envT$res;       
		tunerVal$alg.currentBest <- envT$bst;         
    tunerVal$bfgs = bfgs;
    tunerVal$bfgs$count=nrow(envT$res);
    tunerVal;  
}


