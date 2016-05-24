#*# --------- demo/demo01cpu.r ---------
#*# This demo shows a simple data mining process (phase 1 of TDMR) for the regression task
#*# CPU (from UCI repository, http://archive.ics.uci.edu/ml/datasets/Computer+Hardware).
#*# The data mining process is in main_cpu.r, which calls tdmRegressLoop and tdmRegress
#*# with Random Forest as the prediction model. 

## load package and set working directory
#library(TDMR);
path <- paste(find.package("TDMR"), "demo01cpu",sep="/");
#path <- paste("../inst", "demo01cpu",sep="/");
oldwd <- getwd();
setwd(path);
source("main_cpu.r");                # contains also readCmdCpu()

source("cpu_00.apd",local=TRUE);     # set opts 

result=main_cpu(opts);

## restore old working directory
setwd(oldwd);

