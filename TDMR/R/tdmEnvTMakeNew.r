######################################################################################
# tdmEnvTMakeNew:
#
#' Construct a new environment envT of class \code{\link{TDMenvir}}.
#'
#' Given the general TDMR settings in \code{tdm}, construct an appropriate environment \code{envT}.
#' This is needed as input for \code{\link{tdmBigLoop}}.
#'
#' @param tdm   a list with general settings for TDMR, see \code{\link{tdmDefaultsFill}}  
#'
#' @return Environment \code{envT},  an object of class \code{\link{TDMenvir}},  containing (among others) the elements
#'      \item{\code{runList}}{ \code{=tdm$runList}  }
#'      \item{\code{spotList}}{ \code{=tdm$spotList}  }
#'      \item{\code{tdm}}{ \code{=\link{tdmDefaultsFill}(tdm)}  }
#'      \item{\code{getBst}}{ accessor function(confFile,nExp,theTuner) into \code{envT$bstGrid}   }
#'      \item{\code{getRes}}{ accessor function(confFile,nExp,theTuner) into \code{envT$resGrid}   }
#'      \item{\code{sCList}}{ list of spotConfig-objects, as many as \code{envT$runList} has elements. Each spotConfig object 
#'          sCList[[k]] contains a list \code{opts} as element, which is read from .apd file specified in \code{envT$runList[k]}.  }
#'
#' @seealso   \code{\link{tdmBigLoop}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de}), Patrick Koch
#' @aliases TDMenvir 
#' @export
######################################################################################
tdmEnvTMakeNew <- function(tdm=NULL) {
  ######################################################################################
  # helper fcts for tdmEnvTMakeNew:
  ######################################################################################
  # checkRoiParams: 
  #     Check whether each .roi file in envT$runList contains the same design parameters, otherwise set wP=FALSE and issue a warning
  checkRoiParams <- function(wP,confFile,sC,envT) {
  #firstRoiNames <<- "test"
        if (wP) {
          if (confFile==envT$runList[1]) {
            firstRoiNames <<- rownames(sC$alg.roi);
          } else {
            roiNames = rownames(sC$alg.roi);
            if (length(roiNames)!=length(firstRoiNames)) {
              wP=FALSE;
            } else {
              if (any(roiNames!=firstRoiNames)) wP=FALSE;
            }
          }
          if (wP==FALSE) warning("The design parameters in different ROI files are different --> TDMR sets tdm$withParams to FALSE (no param columns in envT$theFinals)");
        }
        wP;
  }
  # function addSRF adds to opts the element opts$srfFile and reads in case opts$SRF.calc==FALSE from this file the
  # element opts$srf (a list with as many data frames as there are response variables in the task of confFile)
  addSRF <- function(confFile,opts) {
      dir.output <- paste(dirname(opts$dir.output),basename(opts$dir.output),sep="/")  # remove trailing "/", if it exists
      if (!file.exists(dir.output)) {
        success = dir.create(dir.output);
        if (!success) stop(sprintf("Could not create dir.output=%s",dir.output));
      }
      opts$srfFile <- paste(paste(dir.output,confFile,sep="/"),"SRF.Rdata",sep=".")
      if (opts$SRF.calc==FALSE) {
        cat("Loading sorted RF importance from file",opts$srfFile,"...\n")
        if (!file.exists(opts$srfFile)) stop(sprintf("SRF.file=%s does not exist",opts$srfFile));
        srf=NULL;                 # to make package-build-checker happy
        load(file=opts$srfFile);  # load srf
        opts$srf=srf;
      }
      opts;
  }
  ######################################################################################
  # end helper fcts for tdmEnvTMakeNew
  ######################################################################################
  
  envT <- new.env();
  envT$bstGrid <- list();
  envT$resGrid <- list();
  envT$roiGrid <- list();
  envT$sCList <- list();
  envT$theFinals <- NULL;
  envT$tdm <- tdm <- tdmDefaultsFill(tdm);
  envT$runList <- tdm$runList;
  envT$spotList <- tdm$spotList;
  envT$wP <- ifelse(is.null(tdm$withParams), length(tdm$runList)==1, tdm$withParams)
  if (is.null(tdm$runList)) stop("tdm$runList is NULL");
  if (length(unique(tdm$runList)) != length(tdm$runList)) {
    stop(paste("There are duplicates in tdm$runList. Please remove them:\n   ",paste(tdm$runList,collapse=" ")));
#    print(tdm$runList);
#    stop(paste("There are duplicates in tdm$runList. Please remove them."));
  }

  envT <- tdmEnvTAddGetters(envT);

  firstRoiNames <- NULL;    # private storage for checkRoiParams

  #
  # do all necessary file reading **before**  branching into tdmBigLoop 
  # (bigLoopStep in tdmBigLoop.r can be in parallel execution branch, where file access might be not possible)
  #
  k=0;
  for (confFile in envT$runList) {
    k=k+1;
    lastChar=substr(tdm$path,nchar(tdm$path),nchar(tdm$path))
    sepChar=ifelse(lastChar=="/","","/")
    pathConfFile = paste(tdm$path,confFile,sep=sepChar);
    if (!file.exists(pathConfFile)) stop(sprintf("Could not find confFile=%s (current dir=%s)",pathConfFile,getwd()));
    envT$sCList[[k]] <- spotGetOptions(srcPath=tdm$theSpotPath,pathConfFile);
    sC <- envT$sCList[[k]]; print(sC$spot.seed);
    envT$wP <- checkRoiParams(envT$wP,confFile,sC,envT);
    pdFile = sC$io.apdFileName;
  	print(pdFile);
    pdFile = paste(tdm$path,pdFile,sep=sepChar);
    if (!file.exists(pdFile)) stop(sprintf("Could not find pdFile=%s (current dir=%s)",pdFile,getwd()));
    opts <- NULL;     # just to make 'R CMD check' happy  (and in case that after sourcing pdFile 'opts' is not there as expected)
  	source(pdFile,local=TRUE);        # read problem design  (here: all elements of list opts)
		if (is.null(opts))
		  stop(sprintf("%s does not define the required object opts.",pdFile));
    if (class(opts)[1]!="tdmOpts") 
      warning("Object opts is not of class tdmOpts. Consider constructing opts with tdmOptsDefaultsSet().");
    #if (tdm$parallelCPUs>1 & sC$spot.fileMode) {
    #  warning(sprintf("%s: Should not have spot.fileMode==TRUE in parallel execution. spot.fileMode is set to FALSE.",confFile));
    #  envT$sCList[[k]]$spot.fileMode=FALSE;
    #}
    if (!is.null(tdm$TST.trnFrac)) opts$TST.trnFrac=tdm$TST.trnFrac;
    if (!is.null(tdm$TST.valiFrac)) opts$TST.valiFrac=tdm$TST.valiFrac;
    if (!is.null(tdm$READ.target)) opts$READ.target=tdm$READ.target;
    if (!is.null(tdm$oFileMode)) opts$fileMode=tdm$oFileMode;

    opts=tdmOptsDefaultsSet(opts,path=tdm$path);
    opts=addSRF(confFile,opts);      # add opts$srfFile and add opts$srf from file in case opts$SRF.calc==FALSE
    checkOpts(opts);

    if (tdm$parallelCPUs>1 & opts$fileMode==TRUE)
      warning("With tdm$parallelCPUs>1 the setting opts$fileMode=TRUE might be problematic. Consider to set tdm$oFileMode=FALSE");
    
    if (tdm$umode[1]=="TST" & opts$TST.kind=="col" & !(opts$MOD.method %in% c("RF","MC.RF")))
      warning(sprintf("%s: Do you really want tdm$umode=='TST' and opts$TST.kind=='col' when opts$MOD.method is not based on RF? %s",
                      confFile,"You will have no validation data!"));
 	  envT$sCList[[k]]$opts=opts;
  }
  class(envT) <- c("TDMenvir","environment");
  
  envT;
}

######################################################################################
# tdmEnvTReadApd:
#
#' Re-read APD-files (list opts) for an environment envT of class \code{\link{TDMenvir}}.
#'
#' Given an already existing envT with envT$sClist filled, re-read the opts settings from APD files whose filenames 
#' are specified in envT$sClist[[k]]$io.apdfilename.
#'
#' @param envT  an environment of class \code{\link{TDMenvir}}
#' @param tdm   a list with general settings for TDMR (can be envT$tdm)  
#'
#' @return modified environment \code{envT},  an object of class \code{\link{TDMenvir}},  where the items
#'      \item{\code{sCList[[k]]$opts}}{ \code{=tdm$runList}  }
#' are modefied        
#'
#' @seealso   \code{\link{tdmEnvTMakeNew}}
#' @author Wolfgang Konen (\email{wolfgang.konen@@fh-koeln.de}), Patrick Koch
#' @export
######################################################################################
tdmEnvTReadApd <- function(envT,tdm) {
  k=0;
  for (confFile in envT$runList) {
    k=k+1;
    sC <- envT$sCList[[k]]; 
    pdFile = sC$io.apdFileName;
  	print(pdFile);
    pdFile = paste(tdm$path,pdFile,sep="");
    if (!file.exists(pdFile)) stop(sprintf("Could not find pdFile=%s (current dir=%s)",pdFile,getwd()));
    opts <- NULL;     # just to make 'R CMD check' happy  (and in case that after sourcing pdFile 'opts' is not there as expected)
  	source(pdFile,local=TRUE);        # read problem design  (here: all elements of list opts)
		if (is.null(opts))
		  stop(sprintf("%s does not define the required object opts.",pdFile));
    if (class(opts)[1]!="tdmOpts") 
      warning("Object opts is not of class tdmOpts. Consider constructing opts with tdmOptsDefaultsSet().");
    if (!is.null(tdm$TST.trnFrac)) opts$TST.trnFrac=tdm$TST.trnFrac;
    if (!is.null(tdm$TST.valiFrac)) opts$TST.valiFrac=tdm$TST.valiFrac;
    if (!is.null(tdm$READ.target)) opts$READ.target=tdm$READ.target;
    if (!is.null(tdm$oFileMode)) opts$fileMode=tdm$oFileMode;

    opts=tdmOptsDefaultsSet(opts,path=tdm$path);
    checkOpts(opts);

    if (tdm$umode[1]=="TST" & opts$TST.kind=="col" & !(opts$MOD.method %in% c("RF","MC.RF")))
      warning(sprintf("%s: Do you really want tdm$umode=='TST' and opts$TST.kind=='col' when opts$MOD.method is not based on RF? %s",
                      confFile,"You will have no validation data!"));
 	  envT$sCList[[k]]$opts=opts;
  }
  
  envT; 
}

# function checkOpts checks whether there are any new variable names in list opts (e.g. due to misspelling in APD file)
# and if so, issues a NOTE and a warning message.
checkOpts <- function(opts) {
  availNames = c("APPLY_TIME","CLS.CLASSWT","CLS.cutoff","CLS.gainmat","data.title","dir.data","dir.output","dir.Rdata"
                ,"dir.txt","DO.GRAPHICS","DO.POSTPROC","EVALFILE","fct.postproc","fileMode","filename","filesuffix","filetest"
                ,"GD.CLOSE","GD.DEVICE","GD.PNGDIR","GD.RESTART","LOGFILE","logFile","MOD.method","MOD.SEED","ncopies","NRUN","PDFFILE"
                ,"PRE.allNonVali","PRE.knum","PRE.MaxLevel","PRE.PCA","PRE.PCA.npc","PRE.PCA.REPLACE","PRE.SFA","PRE.SFA.doPB"
                ,"PRE.SFA.fctPB","PRE.SFA.npc","PRE.SFA.ODIM","PRE.SFA.PPRANGE","PRE.SFA.REPLACE","PRE.Xpgroup","READ.CMD","READ.INI"
                ,"READ.NROW","READ.TST","READ.TXT","rep","RF.mtry","RF.mtry","RF.nodesize","RF.ntree","RF.OOB","RF.p.all","RF.samp"
                ,"rgain.string","rgain.type","srf","SRF.calc","SRF.cutoff","SRF.kind","SRF.maxS","SRF.method","SRF.minlsi","SRF.ndrop"
                ,"SRF.nkeep","SRF.ntree","SRF.samp","SRF.scale","SRF.verbose","SRF.XPerc","srfFile","SVM.coef0","SVM.cost","SVM.degree","SVM.epsilon"
                ,"SVM.gamma","SVM.kernel","SVM.tolerance","ADA.coeflearn","ADA.mfinal","ADA.rpart.minsplit","test2.show","test2.string"
                ,"TST.COL","TST.kind","TST.NFOLD","TST.SEED","TST.trnFrac","TST.valiFrac","VERBOSE");
  newNames = setdiff(names(opts),availNames);
  if (length(newNames)>0) {
    if (length(newNames)==1) cat("NOTE: A new variable has been defined for list opts: ",newNames,".\n");
    if (length(newNames)>1) cat("NOTE: New variables have been defined for list opts: ",paste(newNames,collapse=", "),".\n");
    warning(paste("NOTE: New variables have been defined for list opts: ",paste(newNames,collapse=", "),"."))
  }
}

######################################################################################
#  tdmEnvTAddGetters:
#'
#' Add getter functions getBst and getRes to environment envT
#'
#' @param   envT       the TDMR environment
#' @return  the augmented \code{envT}
#' @export
######################################################################################
tdmEnvTAddGetters <- function(envT) {  
  # envT$getInd: private helper fct for envT$getBst and envT$getRes
  envT$getInd <- function(envT,confFile,nExp,theTuner) {
    indTuner = which(envT$tdm$tuneMethod==theTuner);
    if (length(indTuner)==0) stop(paste("Could not find tuner ",theTuner,"in envT$tdm$tuneMethod"));
    nConf = which(envT$runList==confFile);
    if (length(nConf)==0) stop(paste("Could not find conf file ",confFile,"in envT$runList"));
    if (nExp<1 | nExp>envT$tdm$nExperim) stop(paste("nExp is not in range {1,...,",envT$tdm$nExperim,"}",sep=""));
    ind = indTuner + length(envT$tdm$tuneMethod)*((nExp-1) + envT$tdm$nExperim*(nConf-1));
  }
  # envT$getBst: return from the linear list envT$bstGrid the data frame bst for the triple {confFile,nExp,theTuner}
  envT$getBst <- function(envT,confFile,nExp,theTuner) {
    ind = envT$getInd(envT,confFile,nExp,theTuner);
    lgrid = length(envT$bstGrid);
    if (ind<1 | ind>lgrid) stop(sprintf("Subscript %d is out of bounds for envT$bstGrid (length is %d)",ind,lgrid));
    envT$bstGrid[[ind]];
  }
  # envT$getRes: return from the linear list envT$resGrid the data frame res for the triple {confFile,nExp,theTuner}
  envT$getRes <- function(envT,confFile,nExp,theTuner) {
    ind = envT$getInd(envT,confFile,nExp,theTuner);
    lgrid = length(envT$resGrid);
    if (ind<1 | ind>lgrid) stop(sprintf("Subscript %d is out of bounds for envT$resGrid (length is %d)",ind,lgrid));
    envT$resGrid[[ind]];
  }
  class(envT) <- c("TDMenvir","environment");
  envT;
}


######################################################################################
#  tdmEnvTLoad:
#' Load an \code{envT}-type environment from file \code{fileRData}. 
#'
#' The loaded envT is augmented with getter functions, see \code{\link{tdmEnvTAddGetters}}.
#'
#' @param   fileRData   string with filename to load. This file is searched in the working dir.
#' @return  envT
#' @export
######################################################################################
tdmEnvTLoad <- function(fileRData) {
  load(fileRData); # loads envT
  envT <- tdmEnvTAddGetters(envT);
  class(envT) <- c("TDMenvir","environment");
  envT;
}

######################################################################################
#  tdmEnvTAddBstRes:
#' Add BST and RES data frames to an existing \code{envT} environment.
#'
#' Load an \code{envT}-type environment from file \code{fileRData}. Its elements 
#' \code{bst}, \code{bstGrid} \code{res}, and \code{resGrid} 
#' overwrite the elements in \code{envT} passed in as argument. 
#'
#' @param   envT       the TDMR environment
#' @param   fileRData   string with filename to load. This file is searched in \code{envT$tdm$path}.
#' @return  the augmented \code{envT}
#' @export
######################################################################################
tdmEnvTAddBstRes <- function(envT,fileRData) {
  loadEnvT <- function(fileRData) {
    load(paste(envT$tdm$path,fileRData,sep="")); # loads envT
    envT;
  }
  envT2 = loadEnvT(fileRData);
  envT$bst = envT2$bst;
  envT$res = envT2$res;
  envT$bstGrid = envT2$bstGrid;
  envT$resGrid = envT2$resGrid;
  if (!is.null(envT2$roiGrid)) envT$roiGrid = envT2$roiGrid;
  if (!is.null(envT2$runList)) envT$runList = envT2$runList;
  if (!is.null(envT2$nExperim)) envT$nExperim = envT2$nExperim;
  if (!is.null(envT2$tuneMethod)) envT$tuneMethod = envT2$tuneMethod;
  class(envT) <- c("TDMenvir","environment");
  envT;
}

######################################################################################
#  tdmEnvTSensi:
#' Make a sensitivity plot based on \code{envT}
#'
#' Given the results from a prior tuning run in \code{envT}, make a sensitivity plot
#' for this run (call SPOT with spotStep="rep")
#' If \code{envT$tdm$nrun > 0} then make additionally with the best-performing parameters from 
#' the tuning run a new unbiased run on the test data.
#'
#' @param   envT  results from a prior tuning run. If envT$tdm$run==0
#' @param   ind   an integer from \code{1:length(envT$bstGrid)}: Take the tuning run with index \code{ind}.
#' @return  \code{finals}, either NULL (no unbiased run) or a one-row data frame with the 
#'          results from the unbiased run
#' @examples            
#'    ## The best results are read from demo02sonar/demoSonar.RData relative to the TDMR 
#'    ## package directory.
#'    oldwd <- getwd(); setwd(paste(find.package("TDMR"), "demo02sonar",sep="/"));
#'    envT = tdmEnvTLoad("demoSonar.RData");    # loads envT
#'    source("main_sonar.r");
#'    envT$tdm$nrun=0;       # =0: don't, >0: do unbiasedRun with opts$NRUN=envT$tdm$nrun
#'    finals = tdmEnvTSensi(envT,1);
#'    if (!is.null(finals)) print(finals);
#'    setwd(oldwd);
#'
#' @export
######################################################################################
tdmEnvTSensi <- function(envT,ind) {
  tdm = envT$tdm;
  if (is.null(tdm$runList)) {
    if (is.null(envT$runList)) stop("Both envT$runList and envT$tdm$runList are NULL");
    tdm$runList = envT$runList;
  }
  nTuner <- length(tdm$tuneMethod);
  nRunList <- length(tdm$runList);
	tuneVec <- rep(tdm$tuneMethod,tdm$nExperim*nRunList);
	expeVec <- rep(sort(rep(1:tdm$nExperim,nTuner)),nRunList);
	confVec <- sort(rep(tdm$runList,tdm$nExperim*nTuner));
  if (ind>length(envT$bstGrid)) stop(sprintf("ind=%d is bigger than length(envT$bstGrid)=%d",ind,length(envT$bstGrid)));
  theTuner = tuneVec[ind];
  nExp = expeVec[ind];
  confFile = confVec[ind];
  nConf = which(confFile==envT$runList);
  #print(c(ind,confFile,nExp,theTuner));
  
  correctEnvTDir(envT,tdm);                # correct direcotories in envT$sCList

  envT$spotConfig <- envT$sCList[[nConf]]; # spotGetOptions(srcPath=tdm$theSpotPath,confFile);  
  envT$bst = envT$bstGrid[[ind]];
  envT$res = envT$resGrid[[ind]];
  #
  # Fix for upgading older BST data frames to the new SPOT version (V1.0.2662): 
  # A BST data frame now has to have a column STEP, otherwise a crash in spotWriteBest (spotHelpFunctions.R) will occur.  
  if (!any(names(envT$bst)=="STEP")) {
    cat("NOTE: Since envT$bst has no column STEP, we add for SPOT a dummy column STEP containing the numbers 1 :",nrow(envT$bst),"\n");
    envT$bst <- cbind(envT$bst,STEP=(1:nrow(envT$bst)));
  }
  envT$spotConfig$alg.roi = envT$roiGrid[[ind]];
  
  envT$spotConfig$spot.fileMode=FALSE;    # call SPOT with spot.fileMode=FALSE (take everything from  envT$spotConfig)
 
  ptm <- proc.time();
  
      # this is for the case spotStep=="rep":
      envT$spotConfig$alg.currentResult <- envT$res;
      envT$spotConfig$alg.currentBest <- envT$bst;	
      
      cat(sprintf("*** Starting SENS REPORT for tuner %s, spotStep=%s, on task %s ***\n",theTuner,"rep",confFile));

      tdmDispatchTuner(theTuner,confFile,"rep",tdm,envT,dataObj);
      # If spotStep=="rep" then tdmDispatchTuner expects that envT$bst, envT$res contain already the right data frames 
      # and it runs just the reporting step of SPOT.
  
  time.txt = c("Proc", "System", "Elapsed");
  if (is.null(tdm$timeMode)) stop("tdm$timeMode is not set (NULL). Consider 'tdm <- tdmDefaultsFill(tdm)' to set all defaults");
  time.TRN=(proc.time()-ptm)[tdm$timeMode]; opts=list(); opts$VERBOSE=1;  
  cat1(opts,paste(time.txt[tdm$timeMode], "time for tuning with tdmDispatchTuner:",time.TRN,"sec\n"));     

  ptm <- proc.time();
  finals <- NULL;
  if (tdm$nrun>0) {
      #
      # read the data and split them into test data and train/vali data
      dataObj <- tdmSplitTestData(envT$spotConfig$opts,tdm,nExp);
      if (!is.null(dataObj)) envT$spotConfig$opts$TST.COL = dataObj$TST.COL;    # this column has to be subtracted in main_* from the input variables

      if (!is.null(tdm$mainFile)) {
        if (!file.exists(tdm$mainFile)) stop(sprintf("Could not find tdm$mainFile=%s (current dir=%s)",tdm$mainFile,getwd()));
 	      source(tdm$mainFile);
      }
    
      for (umode in tdm$umode) {
        cat("*** starting",tdm$unbiasedFunc,"for",confFile,"with umode=",umode,"***\n");
        cmd=paste("finals <-",tdm$unbiasedFunc,"(confFile,envT,dataObj,umode=umode,finals=finals,withParams=envT$wP,tdm=tdm)",sep="");
        eval(parse(text=cmd));
      }
      time.TST=(proc.time()-ptm)[tdm$timeMode];   
      cat1(opts,paste(time.txt[tdm$timeMode], "time for",tdm$unbiasedFunc,":",time.TST,"sec\n"));   
      finals <- cbind(finals,Time.TST=time.TST,Time.TRN=time.TRN);        
      flush.console(); 
  } # if(tdm$nrun>0)
   
  finals;
} # end of function tdmEnvTSensi

######################################################################################
# helper fct for tdmEnvTSensi:
#     correct the path in directories of envT$sCList
correctEnvTDir <- function(envT,tdm) {
    if (!is.null(tdm$path)) {
      for (i in 1:length(envT$sCList)) {    
        # replace the tdm$path-part of opts$dir.* with "./". This workaround necessary when the 
        # original tuning job was done on another machine, e.g. maanbs03, where the path is 
        # expanded to /home/wolfgang/fiwa/... and we run this function on a Windows system
        envT$sCList[[i]]$opts$dir.data = sub(tdm$path,"./",envT$sCList[[i]]$opts$dir.data);
        envT$sCList[[i]]$opts$dir.txt = sub(tdm$path,"./",envT$sCList[[i]]$opts$dir.txt);
        envT$sCList[[i]]$opts$dir.Rdata = sub(tdm$path,"./",envT$sCList[[i]]$opts$dir.Rdata);
        envT$sCList[[i]]$opts$dir.output = sub(tdm$path,"./",envT$sCList[[i]]$opts$dir.output);
      }     
    }
}


