##
## use demo(... , ask=F) to run the demo all at once
##
## Load package
require("SPOT")
## get path of test project
testPath<-find.package("SPOT")
testPath<-file.path(testPath,"demo17Rgp")
## show path
testPath
## Path to demo configuration:
testFile<-file.path(testPath,"rgp0001.conf")
## run example
spotConfig=spot(testFile)

