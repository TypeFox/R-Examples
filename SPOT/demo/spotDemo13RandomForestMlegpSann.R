##
## use demo(... , ask=F) to run the demo all at once
##
## Load package
require("SPOT")
## get path of test project
testPath<-find.package("SPOT")
testPath<-file.path(testPath,"demo13RandomForestMlegpSann")
## show path
testPath
## run example
testFile<-file.path(testPath,"demo13RandomForestMlegpSann.conf")
spotConfig=spot(testFile)
