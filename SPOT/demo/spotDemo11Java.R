##
## use demo(... , ask=F) to run the demo all at once
##
## Load package
require("SPOT")
## get path of test project
testPath<-find.package("SPOT")
testPath<-file.path(testPath,"demo11Java")
## show path
testPath
## This is an example with a user defined custom target function, therefore:
## load the target function interface to workspace
source(file.path(testPath,"bin","spotAlgStartOnePlusOneEsJava.R"))
## run example
testFile<-file.path(testPath,"demo11Java.conf")
spotConfig=spot(testFile)
