##
## use demo(... , ask=F) to run the demo all at once
##
## Load package
require("SPOT")
## get path of test project
testPath<-find.package("SPOT")
testPath<-file.path(testPath,"demo05TgpBranin")
## show path
testPath
## run example
testFile<-file.path(testPath,"demo05TgpBranin.conf")
spotConfig=spot(testFile)
