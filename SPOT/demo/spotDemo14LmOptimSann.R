##
## use demo(... , ask=F) to run the demo all at once
##
## Load package
require("SPOT")
## get path of test project
testPath<-find.package("SPOT")
testPath<-file.path(testPath,"demo14LmOptimSann")
## show path
testPath
## run example
testFile<-file.path(testPath,"demo14LmOptimSann.conf")
spotConfig=spot(testFile)
