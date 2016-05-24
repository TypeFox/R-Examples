######################################################################################
# unbiasedRun
#
#'     Perform unbiased runs with best-solution parameters.
#'
#'     Read the best solution of a parameter-tuning run from \code{envT$bst}, 
#'     execute with these best parameters the function \code{tdm$mainFunc} (usually a classification or 
#'     regression machine learning task), to see 
#'     whether the result quality is reproducible on independent test data or on independently trained models.
#'
#'   @param confFile    the .conf filename, e.g. "appAcid_02.conf"
#'   @param envT        environment, from which we need the objects
#'     \describe{
#'     \item{\code{bst}}{ data frame containing best results (merged over repeats)}
#'     \item{\code{res}}{ data frame containing all results}
#'     \item{\code{theTuner}}{ ["spot"] string}
#'     \item{\code{spotConfig}}{ [NULL] a list with SPOT settings. If NULL, try to read spotConfig from confFile.} 
#'     }
#'   @param dataObj     [NULL] contains the pre-fetched data from which we use here the test-set part. 
#'                      If NULL, set it to \code{tdmSplitTestData(opts,tdm)}  
#'   @param finals      [NULL] a one-row data frame to which new columns with final results are added
#'   @param umode       [ "RSUB" (default) | "CV" | "TST" | "SP_T" ], how to divide in training and test data for the unbiased runs:
#'     \describe{
#'     \item{\code{"RSUB"}}{ random subsampling into (1-tdm$TST.testFrac)\% training and tdm$TST.testFrac\% test data}
#'     \item{\code{"CV"}}{ cross validation (CV) with tdm$nrun folds}
#'     \item{\code{"TST"}}{ all data in opts$filename are used for training, all data in opts$filetest for testing}
#'     \item{\code{"SP_T"}}{ 'split_test': prior to tuning, the data set was split by random subsampling into tdm$TST.testFrac\% test and
#'                  (1-tdm$TST.testFrac)\% training-vali data, tagged via column "tdmSplit". Tuning was done on training-vali data.
#'                  Now we use column "tdmSplit" to select the test data for unbiased evaluation. Training during unbiased evaluation
#'                  is done on a fraction tdm$TST.trnFrac of the training-vali data}
#'     }
#'   @param withParams  [FALSE] if =TRUE, add columns with best parameters to data frame \code{finals}
#'                      (should be FALSE, if different runs have different parameters)
#'   @param tdm         a list with TDM settings from which we use here the elements
#'     \describe{
#'     \item{mainFunc}{ the function to be called for unbiased evaluations}
#'     \item{mainFile}{ change to the directory of mainFile before starting mainFunc}
#'     \item{nrun}{ [5] how often to call the unbiased evaluation}
#'     \item{nfold}{ [10] how many folds in CV (only relevant for umode="CV") }
#'     \item{TST.testFrac}{ [0.2] test set fraction (only relevant for umode="RSUB" or ="SP_T") }
#'     }
#'     The defaults in '[...]' are set by  \code{\link{tdmDefaultsFill}}, if they are not defined on input.
#'   @return \code{finals}     a one-row data frame with final results
#'
#' @note Side Effects:
#'    The list \code{result}, an object of class \code{\link{TDMclassifier}} or \code{\link{TDMregressor}} as returned 
#'    from \code{tdm$mainFunc} is written onto \code{envT$result}. \cr
#'    If  \code{envT$spotConfig} is NULL, it is constructed from confFile. \cr
#'    If \code{spotConfig$opts} (list with all parameter settings for the DM task) is NULL, we try to read it from  
#'    \code{spotConfig$io.apdFileName}. This will issue a warning '... might not work in parallel mode', but is perfectly fine for non-parallel mode. 
#'
#' @examples
#'    ## Load the best results obtained in a prior tuning for the configuration 
#'    ## "sonar_04.conf" with tuning method "spot". The result envT from a prior  
#'    ## run of tdmBigLoop with this .conf is read from demo02sonar/demoSonar.RData.
#'    ## Run task main_sonar again with these best parameters, using the default 
#'    ## settings from tdmDefaultsFill: umode="RSUB", tdm$nrun=5  and tdm$TST.testFrac=0.2.
#'    oldwd <- getwd();          
#'    ## The best results are read from demo02sonar/demoSonar.RData relative to 
#'    ## the TDMR package directory.
#'    setwd(paste(find.package("TDMR"), "demo02sonar",sep="/"));
#'    load("demoSonar.RData");
#'    source("main_sonar.r");
#'    finals <- unbiasedRun("sonar_04.conf",envT,tdm=envT$tdm);
#'    print(finals);
#'    setwd(oldwd);
#'    
#    --- This example is now deprecated, because we do not support .bst-files any longer ----
#    \dontrun{
#    ## If you do not have 'envT' but only a .bst file from a prior tuning run:
#    ## The best results are read from demo02sonar/spot/sonar_04.bst relative to the TDMR package directory.
#    ## (This example is not run automatically, because sonar_04.bst is not in the package distribution) 
#    setwd(paste(find.package("TDMR"), "demo02sonar",sep="/"));
#    envT <- new.env();
#    tdm <- list(mainFunc="main_sonar", tuneMethod="spot");
#    source("main_sonar.r");
#    finals <- unbiasedRun("sonar_04.conf",envT,tdm=tdm);
#    print(finals);
#    setwd(oldwd);
#    }
#'
#' @seealso   \code{\link{tdmBigLoop}}, \code{\link{TDMclassifier}}, \code{\link{TDMregressor}} 
#' @author Wolfgang Konen, FHK, 2010 - 2013
#' @export
#
######################################################################################
unbiasedRun <- function(confFile,envT,dataObj=NULL,finals=NULL,umode="RSUB",withParams=FALSE,tdm=NULL){
    tdm <- tdmDefaultsFill(tdm);
    if (is.null(envT$spotConfig)) envT$spotConfig <- spotGetOptions(srcPath=tdm$theSpotPath,confFile);
    if (is.null(envT$theTuner)) envT$theTuner <- "spot";
    if (is.null(envT$nExp)) envT$nExp <- 1;
    if (is.null(tdm$map)) tdm <- tdmMapDesLoad(tdm); 
    #if (is.null(envT$map)) tdmMapDesLoad(envT,tdm); 
    

    if (is.null(envT$spotConfig$opts)) {
      if (is.null(envT$sCList[[1]]$opts)) {
        pdFile = envT$spotConfig$io.apdFileName;
        warning(paste("List envT$spotConfig does not have the required variable 'opts'."
                     ,"We try to construct it from",pdFile,"(might not work in parallel mode!)")); 
        source(pdFile,local=TRUE);
        opts=tdmOptsDefaultsSet(opts);
      } else {
        opts <- envT$sCList[[1]]$opts;
      }
    } else {
      opts <- envT$spotConfig$opts;
    }
  	opts$ALG.SEED <- envT$spotConfig$alg.seed;

    envT$tdm <- tdm;    # just as information to the caller of this function

    writeLines(paste("start unbiased run for function \"", tdm$mainFunc,"\" ...",sep=""), con=stderr());
# --- this is now in tdmBigLoop.r (before optional parallel execution) ---
#   pdFile = envT$spotConfig$io.apdFileName;
#  	source(pdFile,local=TRUE)         # read problem design  (here: all elements of list opts)   
#   if (!is.null(tdm$mainFile)) source(tdm$mainFile)

  	if (opts$READ.INI) {
  	  if (is.null(dataObj)) dataObj <- tdmSplitTestData(opts,tdm);
      tdm$tstCol=dataObj$TST.COL;  # needed for tdmMapOpts
      #dset <- dataObj$dset;    
      dset <- dsetTrnVa(dataObj);
      tset <- dsetTest(dataObj);       # returns NULL in case of (tdm$umode %in% c("RSUB","CV")
    } else {
      # this branch is deprecated 
      if (is.null(dataObj))
        cat("NOTE: setting opts$READ.INI==FALSE together with argument dataObj==NULL in unbiasedRun.r is deprecated. Consider opts$READ.INI=TRUE.");  
      dset <- NULL;
      tset <- NULL;
    }
    # set certain elements of opts which control selection of training & test data 
    # (depending on switch tdm$umode, see tdmMapDesign.r)
  	tdmOpts <- tdmMapOpts(umode,opts,tdm);           # sets opts$NRUN = tdm$nrun

    # tdmGetObj is deprecated, only for downward compatibility. For actual runs, a simple
    #   bst <- envT$bst
    #   res <- envT$res
    # would be sufficient.
#'     If envT$bst or envT$res is NULL, try to read it from the file (the filename is
#'     inferred envT$spotConfig. If this is NULL, it is constructed from confFile). 
#'     We try to find the files for envT$bst or envT$res in dir envT$theTuner).
    bst <- tdmGetObj(envT$bst,envT$spotConfig$io.bstFileName,envT$theTuner,tdm);
    res <- tdmGetObj(envT$res,envT$spotConfig$io.resFileName,envT$theTuner,tdm);

    k <- nrow(bst);       # last line has the best solution
    #bst <- tdmMapCutoff(bst,k,envT$spotConfig);  # enforce CUTOFF parameter constraint if CUTOFF2[,3,4] appears in .des-file
  	tdmOpts <- tdmMapDesApply(bst,tdmOpts,k,envT$spotConfig,tdm);
  	if (is.null(tdmOpts$fileMode)) tdmOpts$fileMode=TRUE;    # might be necessary for older opts from file

  	cat("Best solution:\n"); print(bst[k,]);
  	#
  	# now tdmOpts has the best solution obtained in a prior tuning, and it has 
  	# training and test set configured according to the current setting of umode

		conf <- as.numeric(bst$CONFIG[k])
		cat(sprintf("Best Config: %5d\n\n",conf))
    #
    # run new model trainings for the configuration stored in tdmOpts (tdm$nrun trainings & evaluations):
    #
		oldwd = getwd();                                               #
    if (!is.null(tdm$mainFile)) setwd(dirname(tdm$mainFile));      # save & change working dir
		result = NULL;        
  	mainCommand <- paste("result <- ", tdm$mainFunc,"(tdmOpts,dset=dset,tset=tset)",sep=" ");
    eval(parse(text=mainCommand));                  # execute the command given in text string mainCommand
    if (!is.list(result)) stop("tdm$mainFunc did not return a list 'result'");
#--- only for testing parallel execution: comment 2 lines above out and uncomment 2 lines below:
#    load("Results2013/DMC2007-base.RData"); 
#    result=envT2$result
#---
		setwd(oldwd);                                                  # restore working dir
    envT$result <- result;

    if (is.null(finals)) {
      # create a data frame with one line of results:
      finals <- data.frame(list(CONF=sub(".conf","",confFile,fixed=TRUE),TUNER=envT$theTuner,NEXP=envT$nExp));
      namFinals <- names(finals);
      if (withParams) {
        pNames=row.names(envT$spotConfig$alg.roi);
        finals <- cbind(finals,bst[k,pNames]);    # bug fix 05/12: this way it works for length(pNames)==1 and for >1    
        names(finals) <- c(namFinals,pNames);     #
      } 
      finals <- cbind(finals
                     , NRUN=tdmOpts$NRUN
                     , NEVAL=nrow(res)
                      );                      
      sgn <- ifelse(class(envT$result)[1]=="TDMregressor",+1,-1);   # minus sign only for classification
      add.finals        <- data.frame(sgn*as.numeric(bst[k,1]),sgn*mean(res$Y));
      names(add.finals) <- paste(tdmOpts$rgain.string,c(".bst",".avg"),sep="");
      finals <- cbind(finals,add.finals);
      #
      # add once the results on the training set from unbiased runs:                      
      suf = ifelse(tdmOpts$MOD.method %in% c("RF","MC.RF"),".OOB",".TRN");
      add.finals <-  data.frame( mean(result$R_train)
                               , sd(result$R_train)
                               );
      names(add.finals) <- paste(c(tdmOpts$rgain.string,"sdR"),suf,sep="");
      finals <- cbind(finals,add.finals);
    } # if(is.null(finals))
      
    # add results on test set from unbiased runs (two columns for each value of umode):
    add.finals <- data.frame(mean(result$R_vali),sd(result$R_vali));
    names(add.finals) <-  paste(c(tdmOpts$rgain.string,"sdR"),umode,sep=".");
    finals <- cbind(finals,add.finals);

    finals;
}

