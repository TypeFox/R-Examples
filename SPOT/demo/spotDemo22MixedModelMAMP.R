##
## use demo(... , ask=F) to run the demo all at once, without interruption
##
## Load package
require("SPOT")
## get path of test project
spotPath<-find.package("SPOT")
demoPath<-file.path(spotPath,"demoMixedModelAnalysis")
## show path
demoPath

## This is a Multiple Algorithm Multiple Problems MAMP example 
## for a mixed model analysis.
## That means, 4 different settings of the Evolution Strategies
## recombination operator are tested on several different 
## instances of the problem. The problem optimized by the ES is 
## the Gaussian Landscape Generator
confFile<-file.path(demoPath,"glges02.conf")
## initialize the design to be evaluated
res<-spot(confFile,spotTask="init")
## run Algorithm and Problem
## please note that this might take several seconds
res<-spot(confFile,spotTask="run")
## generate a report
res<-spot(confFile,spotTask="rep")


