######################################################################################
# Collects garbage
# calls the R garbage collection function gc()
#
collectGarbage <- function()
{
  	while (gc()[2,4] != gc()[2,4]){}
}

######################################################################################
#
# source.tdm: function to load TDMR from source files and to load SPOT and rSFA
#
source.tdm <- function(tdmPath, tdmParallelCPUs=1,theSpotPath=NA,theRsfaPath=NA) {

  createSourcePath <- function(sourceFileName){
    normalizePath(paste(tdmPath,"R",sourceFileName, sep="/"));
  }

  if (!is.na(theSpotPath)) {
      if (theSpotPath=="USE.SOURCE") {
        # load SPOT from source files at a specified location. 
        # (this is useful for debugging SPOT code or for testing new developper versions)
        # --- may need adjustment to your specific SPOT source directory location ---
        if (.Platform$OS.type=="unix") {
          theSpotPath <- '~/Desktop/FH_Koeln/svnspot/trunk/SPOT/R/';
          if (Sys.info()["user"]=="konen") theSpotPath <- '~/svnspot/trunk/SPOT/R/';
        }
        if (.Platform$OS.type=="windows") {
          theSpotPath <- 'C:/WUTemp/FH-MassenDaten/svnspot/trunk/SPOT/R/';
          if (Sys.info()["user"]=="wolfgang")  {
            #theSpotPath <- 'C:/WUTemp/FH-MassenDaten/svnspot/trunk/SPOT/R/';
            theSpotPath <- 'C:/WUTemp/FH-MassenDaten/svnspot/trunk/SPOTsrc/';
          }
        }
      }
      # if 'theSpotPath' is a string different from "USE.SOURCE", then we try 
      # to read the SPOT sources from this path 'theSpotPath'
  }
      
  if (!is.na(theRsfaPath)) {
      if (theRsfaPath=="USE.SOURCE") {
        # load rSFA from source files at a specified location. 
        # (this is useful for debugging rSFA code or for testing new developper versions)
        # --- may need adjustment to your specific rSFA source directory location ---
        if (.Platform$OS.type=="unix") {
          theRsfaPath <- '~/Desktop/FH_Koeln/svnspot/trunk/SPOT/R/';
          if (Sys.info()["user"]=="konen") theRsfaPath <- '~/svnspot/trunk/SPOT/R/';
        }
        if (.Platform$OS.type=="windows") {
          theRsfaPath <- 'C:/WUTemp/FH-MassenDaten/fiwa_soma/trunk/doc/CaseStudies.d/201112.d/tdmExtensions.d/CodeR.d/rSFA/R';
        }
      }
      # if 'theRsfaPath' is a string different from "USE.SOURCE", then we try 
      # to read the rSFA sources from this path 'theRsfaPath'   
  }

  tdmParallel = (tdmParallelCPUs>1);
  if(!tdmParallel){
    require("randomForest");
    require("e1071");        # svm(), Naive Bayes
    #require("matlab");      # repmat() etc., for tdmParaBootstrap.r  - now deprecated 12/2011
    if (is.na(theSpotPath)) {
        require("SPOT");     # load SPOT from the installed library (package version)
    } else {
        cat("Sourcing SPOT from R files in",theSpotPath,"\n");    
        sourceSpotFromDir <- function(sdir) {
            oldwd=getwd(); setwd(sdir);
            for (f in dir())   {
              if (file.info(f)$isdir) {
                sourceSpotFromDir(f);
              } else {            
                source(f);
              }
            }
            setwd(oldwd);
        }
        sourceSpotFromDir(theSpotPath);         # recursive call (march through all subdirs)
    }
    
    if (is.na(theRsfaPath)) {
        require("rSFA");     # load rSFA from the installed library (package version)
    } else {
        oldwd=getwd(); setwd(theRsfaPath);
        for (f in dir())   source(f);
        setwd(oldwd);
    }
    
    source(createSourcePath("makeTdmRandomSeed.r"))
    source(createSourcePath("printTDMclassifier.r"))
    source(createSourcePath("printTDMregressor.r"))
    source(createSourcePath("tdmClassify.r"))
    source(createSourcePath("tdmClassifyLoop.r"))
    source(createSourcePath("tdmEmbedDataFrame.r"))
    source(createSourcePath("tdmGeneralUtils.r"))
    source(createSourcePath("tdmGraphicUtils.r"))
    source(createSourcePath("tdmMetacostRf.r"))
    source(createSourcePath("tdmModelingUtils.r"))    
    source(createSourcePath("tdmOptsDefaults.r"))
    source(createSourcePath("tdmParaBootstrap.r"))
    source(createSourcePath("tdmPreprocUtils.r"))
    source(createSourcePath("tdmReadData.r"))
    source(createSourcePath("tdmRegress.r"))
    source(createSourcePath("tdmRegressLoop.r"))

    source(createSourcePath("tdmBigLoop.r"))
    source(createSourcePath("tdmDefaultsFill.r"))
    source(createSourcePath("tdmDispatchTuner.r"))
    source(createSourcePath("tdmEnvTMakeNew.r"))
    source(createSourcePath("tdmExecSpotStep.r"))
    source(createSourcePath("tdmGetObj.r"))
    source(createSourcePath("tdmMapDesign.r"))
    source(createSourcePath("tdmPlotResMeta.r"))
    source(createSourcePath("tdmROCR.r"))
    source(createSourcePath("tdmSplitTestData.r"))
    source(createSourcePath("tdmStartSpot.r"))
    source(createSourcePath("tdmStartOther.r"))
    source(createSourcePath("unbiasedRun.r"))
    source(createSourcePath("unbiasedBestRun_O.r"))
    
  }
  else   # i.e. if (tdmParallel)
  {
    stop("For tdm$parallelCPUs>1 it is required to use the *library* version of TDMR. Consider to set tdm$tdmPath=NULL.")
    #sourceTDMR_SF(theSpotPath,theRsfaPath)        # --- obsolete now, the old source is in R-DM-Template-deprecated\..\prepareParallelExec.r
  }
    
  collectGarbage()
}


