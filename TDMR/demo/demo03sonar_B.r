#*# --------- demo/demo03sonar_B.r ---------
#*# Same as demo/demo03sonar.r, but with parameters for multiple tuning experiments & longer tuning runs:
#*#    in demo03sonar_B.r:    tdm$nExperim=2; tdm$nrun=5;
#*#    and in sonar_05.conf:  auto.loop.nevals = 50; init.design.size = 10;
#*# and with two tuners SPOT and LHD in comparison.

## load package and set working directory (dir with .apd, .conf and main_*.r file)
path <- paste(find.package("TDMR"), "demo02sonar",sep="/");
#path <- paste("../inst", "demo02sonar",sep="/");

## preliminary settings for TDMR
tdm <- list( mainFunc="main_sonar"
            , runList = "sonar_05.conf"
            , umode=c("RSUB")           # ["CV" | "RSUB" | "TST" | "SP_T" ]
            , tuneMethod = c("spot","lhd")
            , filenameEnvT="demo03.RData"   # file to save environment envT (in working dir)
            , nrun=5, nfold=2         # repeats and CV-folds for the unbiased runs
            , nExperim=2
            , parallelCPUs=1
            , parallelFuncs="readCmdSonar"  # fct's to export in addition to tdm$mainFunc in case parallelCPUs>1
            , optsVerbosity = 0       # the verbosity for the unbiased runs
            );
## tdm$runList="sonar_05.conf" has the settings for the tuning process (e.g. 
##    - "auto.loop.steps"=number of SPOT generations       
##    - "auto.loop.evals"=budget of model building runs and 
##    - io.roiFileName = "sonar_05.roi"
## ). tdm$runList could contain other files as well (e.g. 
##    c("sonar_01.conf","sonar_02.conf","sonar_03.conf")
## ), if desired.

spotStep = "auto";    ## spotStep can be either "auto" (do automatic tuning) or 
                      ## "rep" (make a visual report and an unbiased run on best results)
source(paste(path,tdm$mainFile,sep="/"));    
source(paste(path,"start_bigLoop.r",sep="/"),chdir=TRUE);    # change dir to 'path' while sourcing

## the resulting tuning surface (the metamodel) can be inspected interactively with
##      tdmEnvTLoad(paste(path,tdm$filenameEnvT,sep="/"));     
##      tdmPlotResMeta(envT);
## (tdmEnvTLoad(...) is only needed for reloading envT in another R-session)
