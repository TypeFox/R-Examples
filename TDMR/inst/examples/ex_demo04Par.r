#*# This demo is the same as demo04cpu.r, except that it is executed on 2 cluster nodes in parallel.
#*# The only differences to demo04cpu.r are:
#*#     a) tdm$parallelCPUs=2
#*#     b) tdm$parallelFuncs=c("readCmdCpu","cpu.postproc")
#*# The latter is only necessary if main_cpu.r contains extra functions besides tdm$mainFunc="main_cpu",
#*# in this case "cpu.postproc" and "readCmdCpu", which need to be clusterExport'ed to the cluster nodes

## load package and set working directory
require(TDMR);
#path <- paste(find.package("TDMR"), "demo01cpu/",sep="/");
path <- paste("../../inst", "demo01cpu/",sep="/");
oldwd <- getwd();  setwd(path);
#
source(paste(path,"main_cpu.r",sep=""));

## preliminary settings for TDMR
tdm <- list(    mainFunc="main_cpu"
#              , path=path
              , umode=c("RSUB")    # ["CV" | "RSUB" | "TST"]
              , tuneMethod=c("spot","lhd")
              , filenameEnvT="demo04cpu.RData"   # file to save environment envT (in working dir)
              , finalFile="cpu.fin"     # where to write final results (best solution & unbiased eval for each tuner/.conf-combination)
              , withParams=TRUE         # list the columns with tuned parameter in final results
              , optsVerbosity=0         # the verbosity for the unbiased runs
              , parallelCPUs=2
              , parallelFuncs=c("readCmdCpu","cpu.postproc")
              );
## fill in other defaults for list tdm
tdm <- tdmDefaultsFill(tdm);
## cpu_01.conf has the settings for the tuning process (e.g. "auto.loop.steps"=number of SPOT generations
## "auto.loop.evals"=budget of model building runs and io.roiFileName = "cpu_01.roi").
tdm$runList = c("cpu_01.conf");
## spotStep can be either "auto" (do automatic tuning) or "rep" (make a visual report and an unbiased run on best results)
spotStep = "auto";

## construct an initial environment envT from the given TDMR settings in tdm
## (this contains also the fill-in of other defaults for tdm via
##      envT$tdm <- tdmDefaultsFill(tdm);
## )
envT <- tdmEnvTMakeNew(tdm);

## the call to tdmBigLoop will start the whole TDMR process:
## - for each file in tdm$runList a complete DM tuning is started with each tuning method tdm$tuneMethod  (if spotStep=="auto")
## - the best result from tuning is fed into an unbiased model build and evaluation run
## - results are printed and returned in envT$theFinals
## - more detailed results are in other elements of environment envT
## - two plots:
##      a) the progression of the response variable Y and the parameter variables during tuning
##      b) the sensitivity plot for each parameter in the vicinity of the best solution found
envT <- tdmBigLoop(envT,spotStep);

## restore old working directory
setwd(oldwd);

## the resulting tuning surface (the metamodel) can be inspected interactively with
##      load(paste(path,tdm$filenameEnvT,sep="/"));
##      tdmPlotResMeta(envT);
## (load(...) is only needed for reloading envT in another R-session)
