##
## use demo(... , ask=F) to run the demo all at once
##
## Load package
require("SPOT")
## get path of test project
testPath<-find.package("SPOT")
testPath<-file.path(testPath,"demo07RandomForestSann")
## show path
testPath
## run example
testFile<-file.path(testPath,"demo07RandomForestSann.conf")
spotConfig=spot(testFile)
