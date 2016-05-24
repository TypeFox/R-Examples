#*# --------- demo/demo08parallel.r ---------
#*# This demo does the same as demo03sonar.r, but it runs 4 experiment on 4 parallel cores
#*# (if your environment supports parallel clusters with the R-core package 'parallel'). 
#*#
#*# Note on your task manager, how for a short time 4 processes are active (usually 'RScript').
#*# 
#*# Tuning runs are rather short, to make the example run quickly. 
#*# Do not expect good numeric results. 
#*# See demo/demo03sonar_B.r for a somewhat longer tuning run, with two tuners SPOT and LHD.

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
            , tuneMethod =  c("spot","lhd") # c("lhd")   #
            , filenameEnvT="demo03.RData"   # file to save environment envT (in working dir)
            , nrun=2, nfold=2         # repeats and CV-folds for the unbiased runs
            , nExperim=2
            , parallelCPUs=4
            , parallelFuncs=c("readCmdSonar")
            , optsVerbosity = 0       # the verbosity for the unbiased runs
            , oFileMode = FALSE       # set opts$fileMode==FALSE
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

## construct an initial environment envT from tdm
## (this contains also tdmDefaultsFill(tdm))
## then run tdmBigLoop
envT <- tdmExecSpotStep(tdm,spotStep);

setwd(oldwd);               # restore old working directory

## the resulting tuning surface(s) (the metamodel(s)) can be inspected interactively with
##      load(paste(path,tdm$filenameEnvT,sep="/"));     
##      tdmPlotResMeta(envT);
## (load(...) is only needed for reloading envT in another R-session)
