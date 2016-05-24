######################################################################################
# tdmMapDesLoad
#
#'    Load the mapping files. 
#'
#'    Load the map files \code{"tdmMapDesign.csv"} and optionally 
#'    also \code{"userMapDesign.csv"} and store them in \code{tdm$map} and \code{tdm$mapUser},
#'    resp. These maps are used by \code{\link{tdmMapDesApply}}.   \cr
#'    \code{"tdmMapDesign.csv"} is searched in the TDMR library path \code{find.package("TDMR")}.
#'    (For the developer version: \code{<tdm$tdmPath>/inst}). \cr
#'    \code{"userMapDesign.csv"} is searched in tdm$path (which is getwd() if the user did not define tdm$path).
#'
#' @param tdm   list, needed for \code{tdm$tdmPath} and \code{tdm$path}
#' @return \code{tdm}, the modified list with new elements \code{tdm$map} and \code{tdm$mapUser}
#' @seealso  \code{\link{tdmMapDesApply}}
#' @export
tdmMapDesLoad <- function(tdm=list()) {
      if (is.null(tdm$tdmPath)) {
        mapPath <- find.package("TDMR");               # this is for the package version
      } else {
        mapPath <- paste(tdm$tdmPath,"inst",sep="/");   # this is for the developer source version
      }
      tdmMapFile=paste(mapPath,"tdmMapDesign.csv",sep="/")
      if (!file.exists(tdmMapFile)) stop(sprintf("Could not find map file %s",tdmMapFile));
      tdm$map <- read.table(tdmMapFile,sep=";",header=T);
      cat(sprintf("Read tdmMapDesign.csv from %s\n",mapPath));
      #userMapDir=ifelse(is.null(tdm$mainFile),".",dirname(tdm$mainFile));
      userMapFile=paste(tdm$path,"userMapDesign.csv",sep="");
      if (file.exists(userMapFile)) {
        tdm$mapUser <- read.table(userMapFile,sep=";",header=T);
        cat(sprintf("Read userMapDesign.csv from %s\n",tdm$path));
      } else {
        tdm$mapUser <- NULL;
        cat(sprintf("Note: No file userMapDesign.csv in %s\n",tdm$path));
      }
      tdm;
}    

# tdmMapDesApply
#
#'    Apply the mapping from \code{des} to \code{opts}.     
#'
#'    For each variable which appears in .roi (and thus in .des file and design point data frame \code{des}): 
#'    set its counterpart in list \code{opts} to the values of the \code{k}-th row in \code{des}.
#'    For each variable not appearing: leave its counterpart in \code{opts} at its default value from .apd file.
#'
#' @param des   design points data frame
#' @param opts  list of options
#' @param k     apply mapping for the \code{k}-th design point
#' @param spotConfig  list, we needed here \code{spotConfig$alg.roi} and \code{envT$mapUser}, see \code{\link{tdmMapDesLoad}},
#'              and in addition \code{envT$spotConfig$alg.roi} 
#' @param tdm   list, we need here \code{tdm$map} and \code{tdm$mapUser}
#' @return \code{opts}, the modified list of options
#' @seealso  \code{\link{tdmMapDesLoad}}
#' @export
tdmMapDesApply <- function(des,opts,k,spotConfig,tdm) {
    cRound <- function(n,map,x) {
      x <- ifelse(map$isInt[n]==1,round(x),x);
      x; 
    }
    setMapValue <- function(n,map,des,opts,k) {
    	cmd = paste("if(!is.null(des$",map$roiValue[n],")) ",map$optsValue[n],"=cRound(n,map,des$",map$roiValue[n],"[k])",sep="");
      # (The check "if (!is.null...)" assures that tdmMapDes will work for all different .roi files.)
    	eval(parse(text=cmd));
    	opts;
    }
    
    # check whether each parameter column in des can be found in map or mapUser (thus whether it can be mapped to a variable in opts):
    pNames=row.names(spotConfig$alg.roi);
    for (d in pNames) {
      if (length(which(tdm$map$roiValue==d))+length(which(tdm$mapUser$roiValue==d))==0)
        stop(sprintf("tdmMapDesApply: cannot find a mapping for design variable %s.\n  ",d),
             sprintf("Please check spelling in ROI file %s %s", spotConfig$io.roiFileName,
                     "or extend tdmMapDesign.csv or userMapDesign.csv appropriately!"));
    }
    
    if (nrow(tdm$map)>0) {
      for (n in 1:nrow(tdm$map)) {
        opts <- setMapValue(n,tdm$map,des,opts,k);
      }
    }
    if (!is.null(tdm$mapUser)) {
      if (nrow(tdm$mapUser)>0) {
        for (n in 1:nrow(tdm$mapUser)) {
          opts <- setMapValue(n,tdm$mapUser,des,opts,k);
        }
      }
    }
    
    # TODO: check and extend this for arbitrary class length
    if (!is.null(opts$CLS.CLASSWT) && is.na(opts$CLS.CLASSWT[1])) opts$CLS.CLASSWT[1]=10;
    
    opts;
}

######################################################################################
tdmMapDesInt <- function(des,printSummary=T,spotConfig=NULL)
{
	# the following is a fix for the fact that the current SPOT version does not
	# generate INTs but FLOATs for each INT in the .roi file.
	# We hard-code a round() for certain names we know that they only can be INT:
	if (!is.null(des$MTRY)) des$MTRY <- round(des$MTRY);
	if (!is.null(des$NTREE)) des$NTREE <- round(des$NTREE);
	if (!is.null(des$NODESIZE)) des$NODESIZE <- round(des$NODESIZE);
	if (!is.null(des$SAMPSIZE1)) des$SAMPSIZE1 <- round(des$SAMPSIZE1);
	if (!is.null(des$SAMPSIZE2)) des$SAMPSIZE2 <- round(des$SAMPSIZE2);
	if (!is.null(des$SAMPSIZE3)) des$SAMPSIZE3 <- round(des$SAMPSIZE3);
	if (!is.null(des$SAMPSIZE4)) des$SAMPSIZE4 <- round(des$SAMPSIZE4);
	if (!is.null(des$SAMPSIZE5)) des$SAMPSIZE5 <- round(des$SAMPSIZE5);
	if (!is.null(des$PRE_NPC)) des$PRE_NPC = round(des$PRE_NPC);
  #
  # --- TODO: infer the INT columns automatic from spotConfig's ROI
  # --- (when the interface to tdmStart* is clear everywhere)

	if (printSummary) print(summary(des));
	des;
}

######################################################################################
# tdmMapOpts:
#   helper fct for unbiasedRun_*.r, map certain parameters of opts for the unbiased
#   runs, depending on parameter umode, which controls the mode of the unbiased run
# INPUT
#   umode       # one of { "RSUB" | "CV" | "TST" | "SP_T" }
#               # ="RSUB": activate random subsampling with tdm$TST.testFrac test data
#               # ="CV": activate cross validation with tdm$nfold [10] folds
#               # ="TST": activate test on unseen test data (user-def'd)
#               # ="SP_T": activate test on unseen test data (randomly selected by tdmSplitTestData)
#   opts        # current state of parameter settings
#   tdm         # list, here we use the elements
#     nfold         # [10] value for opts$TST.NFOLD during unbiased runs with umode="CV"
#     tstCol        # ["TST"] value for opts$TST.COL during unbiased runs with umode="TST"
#     TST.testFrac  # [0.2] value for opts$TST.valiFrac during unbiased runs with umode="RSUB"
#     U.trnFrac     # [NULL] value for opts$TST.trnFrac during unbiased runs with umode="SP_T"  or umode="TST"
#     nrun          # [5] value for opts$NRUN during unbiased runs
#     test2.string  # ["default.cutoff"] value for opts$test2.string during unbiased runs
#   (the defaults in '[...]' are filled via tdmDefaultsFill, if tdm does not define them.)
# OUTPUT
#   opts        #  enhanced state of parameter settings
######################################################################################
tdmMapOpts <- function(umode,opts,tdm)
{
    tdm <- tdmDefaultsFill(tdm);      # fill in default values if they are not yet set
    
    setOpts.RSUB <- function(opts) {
      opts$TST.kind <- "rand" # select test set by random subsampling, see tdmModCreateCVindex in tdmModelingUtils.r
      opts$TST.valiFrac=tdm$TST.testFrac # set this fraction of data aside for testing 
      opts;
    } 
    setOpts.TST <- function(opts) {
      opts$TST.kind <- "col"  # select test set from column TST.COL, see tdmModCreateCVindex in tdmModelingUtils.r
      opts$READ.TST = T       # =T: read in extra unseen test data 
                              # and fill column dset[,opts$TST.COL] accordingly (1 for test records) 
      if (!is.null(tdm$U.trnFrac)) 
        opts$TST.trnFrac=tdm$U.trnFrac  # Use only this fraction of the TrainVa-data for training (even when opts$TST.kind="col").
                                        # If tdm$U.trnFrac is NULL, then tdmModCreateCVindex uses all TrainVa-data for training.
      opts$TST.COL=tdm$tstCol                              
      opts;
    } 
    setOpts.CV <- function(opts) {
      opts$TST.kind <- "cv"       # select test data by CV, see tdmModCreateCVindex in tdmModelingUtils.r
      opts$TST.NFOLD = tdm$nfold  # number of CV-folds 
      opts$READ.TST = F           # =F: do not read in extra unseen test data 
      opts;
    } 
    setOpts.SP_T <- function(opts) {
      # [we leave opts$TST.kind at its current value]
      if (!is.null(tdm$U.trnFrac)) 
        opts$TST.trnFrac=tdm$U.trnFrac  # Set this fraction of the TrainVa-data aside for training (only for opts$TST.kind="rand").
                                        # If tdm$U.trnFrac is NULL, then tdmModCreateCVindex sets opts$TST.trnFrac = 1-opts$TST.valiFrac.
      if (!is.null(tdm$TST.valiFrac)) 
        opts$TST.valiFrac=tdm$TST.valiFrac
      opts$TST.NFOLD = tdm$nfold  # number of CV-folds (only for opts$TST.kind=="cv")
      opts;
    } 
    opts <- switch(umode
      , "RSUB" = setOpts.RSUB(opts)
      , "CV" = setOpts.CV(opts)
      , "TST" = setOpts.TST(opts)
      , "SP_T" = setOpts.SP_T(opts)
      , "INVALID"
      );
    if (opts[1]=="INVALID") 
      stop(sprintf("*** Invalid umode=%s ***\n",umode));

    opts$NRUN = tdm$nrun    # how many runs with different train & test samples  - or -
                            # how many CV-runs, if TST.kind=="cv"
    opts$test2.string=tdm$test2.string   # test2, another copy of the test set, receives a special treatment and the
                            # special treatment is either opts$test2.string = "no postproc" or = "default cutoff"
#    opts$READ.TXT = T       # =T: read data from .csv and save as .Rdata, =F: read from .Rdata

    opts$VERBOSE <- opts$SRF.verbose <- tdm$optsVerbosity;
    
    opts;
}

#--- DEPRECATED, just for compatibility with older versions unbiasedBestRun_appStorm.R, unbiasedBestRun_DMC2007.R -----------------
#
tdmMapOpts.OLD <- function(umode,opts,test2.string="no postproc",nfold=10) {
  tdm = list();
  tdm$test2.string=test2.string;
  tdm$nfold=nfold;
  opts <- tdmMapOpts(umode,opts,tdm);
}

######################################################################################
# ---- tdmMapCutoff is deprecated, all is done in tdmModAdjustCutoff now. ---
# tdmMapCutoff:
#   helper fct for tdmStart*.r in the classification task case
#   enforce parameter constraint sum(CUTOFFi)=1 if CUTOFF2,3,4 appears in .des-file.
#   The following code works up to N.c=5  [N.c = number of classes]
#
#   Bug fix 01/2011: the constraint fix is refined in such a way that all constrained will 
#   be enforced simultaneously: 
#       1) sum(cutoffs)==1
#       2) any(cutoffs) <= high ROI constraint (from spotConfig)
#       3) any(cutoffs) >= low ROI constraint (from spotConfig)
#   We assume here for simplicity that the high ROI and low ROI are the same for all cutoffs.
######################################################################################
tdmMapCutoff <- function(des,k,spotConfig) {   
    # note that nothing needs to be done here if only des$CUTOFF1 appears (two-class problem):
    # then the dependent parameter opts$CLS.cutoff[2] = 1-opts$CLS.cutoff[1] is set in tdmModAdjustCutoff
    #
  	if (!is.null(des$CUTOFF2)) {        # enforce parameter constraint if CUTOFF1,2,3,4 appears in .des-file:
      hig=spotConfig$alg.roi["CUTOFF1","upper"];
      low=spotConfig$alg.roi["CUTOFF1","lower"];    
      OLD.VER=F;	 # if OLD.VER==T: restore the previous version which had the low/high ROI constraint
			             # violation bug (but gave better results on appAcid)	 
      constr=ifelse(OLD.VER,0,1);
  	  # If sum(cutoff[1:(N.c-1)])>=1, we clip cutoff[1:(N.c-1] in such a way that the sum is below 1. 
      # Why eps? - Because RF does not allow cutoff[N.c]=0, it has to be positive (see main_DMC2007.r)
      eps <- runif(1, min=0.0001, max=0.03);
      csum <- des$CUTOFF1[k]+des$CUTOFF2[k];
      if (!is.null(des$CUTOFF3)) csum <- csum + des$CUTOFF3[k];
      if (!is.null(des$CUTOFF4)) csum <- csum + des$CUTOFF4[k];
      if (csum+eps>constr) {
        # If constr==1: CUTOFF1..4 sum up to something larger than 1-eps.
        # There is no room left for CUTOFF[N.c], so we scale CUTOFF1..4 down.
        # If constr==0 (i.e. OLD.VER==T): this if-branch is always true, we scale always. 
        # This means that CUTOFF[N.c]=eps and that CUTOFF1..4 will be scaled up, if csum+eps<1; 
        # in this case the high ROI constraint can be violated.        
        des$CUTOFF1[k] <- des$CUTOFF1[k]/(csum+eps);
        des$CUTOFF2[k] <- des$CUTOFF2[k]/(csum+eps);
        if (!is.null(des$CUTOFF3)) des$CUTOFF3[k] <- des$CUTOFF3[k]/(csum+eps);
        if (!is.null(des$CUTOFF4)) des$CUTOFF4[k] <- des$CUTOFF4[k]/(csum+eps);
      }
      # CUTOFF[N.c] will be set in main_TASK.r
      # Attention: it can (and will) happen, that 1-csum>hig, if hig is the high ROI 
      # constraint for CUTOFF. In this case CUTOFF[N.c]=1-csum will be larger than hig.
      # We currently tolerate this; CUTOFF[N.c]=1-csum will still obey the constraint "<1". 
      
      cutoffs=c(des$CUTOFF1[k],des$CUTOFF2[k],des$CUTOFF3[k],des$CUTOFF4[k]);
      # cutoffs will be shorter than 4 elements if some of the des$CUTOFF*[k] are NULL
      
      if (!OLD.VER) {
        while (any(cutoffs<low)) {
          # repair violations of the low ROI constraint, if any:
          # add the missing amount delta to the lowest cutoff, subtract it from the highest
          cutoffs = c(cutoffs,1-sum(cutoffs));
          wmin = which(cutoffs==min(cutoffs));
          wmax = which(cutoffs==max(cutoffs));
          delta = low-cutoffs[wmin];
          cutoffs[wmin] = cutoffs[wmin]+delta;
          cutoffs[wmax] = cutoffs[wmax]-delta;
          cutoffs=cutoffs[1:(length(cutoffs)-1)];
        }
        if (!is.null(des$CUTOFF1[k])) des$CUTOFF1[k]=cutoffs[1];
        if (!is.null(des$CUTOFF2[k])) des$CUTOFF2[k]=cutoffs[2];
        if (!is.null(des$CUTOFF3[k])) des$CUTOFF3[k]=cutoffs[3];
        if (!is.null(des$CUTOFF4[k])) des$CUTOFF4[k]=cutoffs[4];
      }
      
      if (any(cutoffs>hig)) {
        msg=sprintf("Some cutoffs violate high ROI constraint %f (for k=%d)",hig,k);
        warning(msg);
        print(cutoffs);
      }
      if (any(cutoffs<low)) {
        msg=sprintf("Some cutoffs violate low ROI constraint %f (for k=%d)",low,k);
        warning(msg);
        print(cutoffs);
      }
    } # if (!is.null(des$CUTOFF2))
#print(c(cutoffs,1-sum(cutoffs)));
#browser()   
    des;
}
