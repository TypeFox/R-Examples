#!/usr/bin/Rscript
# vim:set ff=unix expandtab ts=2 sw=2:
# extend Runit to check warnings
library("RUnit")
library("deSolve")
library("parallel")
source("testhelpers.R")
dataPrefix="../../data/"
dataPaths=Sys.glob(paste(dataPrefix,"*.rda",sep=""))
for (dfile in dataPaths){
  load(dfile)
}
#load("../../data/C14Atm_NH.rda")
#load("../../data/HarvardForest14CO2.rda")
#load("../../data/IntCal09.rda")
#load("../../data/afn_Atm_NH.rda")
prefix="../../R/"
globstring=paste(prefix,"*.R",sep="")
auto_paths=Sys.glob(globstring)
preload_paths=c(
         #  "1_TimeMap.R",
          # "Model.R"
          # "Model_14.R",
          # "DecompositionOperator.R",
          # "GeneralModel14.R",
          )
preload_paths=sapply(preload_paths,function(x){return(paste(prefix,x,sep=""))})
print(preload_paths)
all_paths=c(preload_paths,auto_paths)

for (f in all_paths){
    print(f)
    source(f,echo=FALSE)
}
