#*# --------- demo/demo03sonar_A.r ---------
#*# This demo is in function identical to demo/demo03sonar.r, it just shows 
#*# the code contents (and the comments) of tdmExecSpotStep explicitely.

## load package and set working directory (dir with .apd, .conf and main_*.r file)
#library(TDMR);
path <- paste(find.package("TDMR"), "demo02sonar",sep="/");
#path <- paste("../inst", "demo02sonar",sep="/");
oldwd <- getwd();  setwd(path);
source("main_sonar.r");    # from working dir

## preliminary settings for TDMR
tdm <- list( mainFunc="main_sonar"
            , runList = "sonar_04.conf"
            , umode=c("CV")           # ["CV" | "RSUB" | "TST" | "SP_T" ]
            , tuneMethod = c("lhd")   # c("spot","lhd")
            , filenameEnvT="demo03.RData"   # file to save environment envT (in working dir)
            , nrun=1, nfold=2         # repeats and CV-folds for the unbiased runs
            , nExperim=1
            , parallelCPUs=1
            , parallelFuncs=c("readCmdSonar")
            , optsVerbosity = 0       # the verbosity for the unbiased runs
            );
## tdm$runList="sonar_04.conf" has the settings for the tuning process (e.g. 
##    - auto.loop.steps = number of SPOT generations       
##    - auto.loop.evals = budget of model building runs and 
##    - io.roiFileName = "sonar_04.roi"
## ). tdm$runList could contain other files as well (e.g. 
##    c("sonar_01.conf","sonar_02.conf","sonar_03.conf")
## ), if desired.

spotStep = "auto";    ## spotStep can be either "auto" (do automatic tuning) or 
                      ## "rep" (make a visual report and an unbiased run on best results)

    ## the contents of tdmExecSpotStep is written out explicitely (just to show the code)
    if (spotStep == "auto") {
      #
      # perform a complete tuning + unbiased eval
      #
      
      ## construct an initial environment envT from the given TDMR settings in tdm
      ## (this contains also the fill-in of other defaults for tdm via
      ##      envT$tdm <- tdmDefaultsFill(tdm);
      ## )
      envT <- tdmEnvTMakeNew(tdm);

      ## the call to tdmBigLoop will start the whole TDMR process:
      ## - for each file in tdm$runList a complete DM tuning is started with each tuning
      ##   method tdm$tuneMethod  (if spotStep=="auto")
      ## - the best result from tuning is fed into an unbiased model build and evaluation run
      ## - results are printed and returned in envT$theFinals
      ## - more detailed results are in other elements of environment envT
      ## - two plots:
      ##      a) the progression of the response variable Y and the parameter variables during tuning
      ##      b) the sensitivity plot for each parameter in the vicinity of the best solution found
      envT <- tdmBigLoop(envT,spotStep);
    } else        
    {        # i.e. spotStep == "rep" or == "report"
      #
      # re-use prior tuning result from tdm$filenameEnvT:
      # do only spot report and unbiased eval on best tuning solution
      #
      if (is.null(tdm$filenameEnvT)) tdm$filenameEnvT=sub(".conf",".RData",tdm$runList[1],fixed=TRUE);
      envT <- tdmEnvTLoad(tdm$filenameEnvT);
      envT <- tdmBigLoop(envT,"rep");
    }

setwd(oldwd);               # restore old working directory

## the resulting tuning surface(s) (the metamodel(s)) can be inspected interactively with
##      load(paste(path,tdm$filenameEnvT,sep="/"));     
##      tdmPlotResMeta(envT);
## (load(...) is only needed for reloading envT in another R-session)
