################################################################################
###
### merge and thin out results on rdafile
###
require(RobAStRDA)
.saveGridToRda <- RobAStRDA:::.saveGridToRda
.readGridFromCSV <- RobAStRDA:::.readGridFromCSV
.mergeGrid <- RobAStRDA:::.mergeGrid
.MakeSmoothGridList <- RobAStRDA:::.MakeSmoothGridList
.computeInterpolators <- RobAStRDA:::.computeInterpolators
.mergeF <- RobAStRDA:::.computeInterpolators
.generateInterpolators <- RobAStRDA:::.generateInterpolators
#####

oldwd <- getwd()
.basepath <- "C:/rtest/RobASt/branches/robast-1.0/pkg"
.myFolderFrom <- file.path(.basepath,"RobExtremesBuffer")
myRDA0 <- file.path(.basepath,"RobExtremesBuffer/sysdata.rda")
#myRDA <- file.path(.basepath,"RobExtremesBuffer/sysdata.rda")
#myRDA0 <- file.path(.basepath,"RobAStRDA/R/sysdata0.rda")
myRDA <- file.path(.basepath,"RobAStRDA/R/sysdata.rda")
CSVFiles <- grep("\\.csv$", dir(.myFolderFrom), value=TRUE)
CSVFiles <- paste(.myFolderFrom, CSVFiles, sep="/")

.saveGridToRda(CSVFiles, toFileRDA = myRDA0, withMerge = FALSE,
               withPrint = TRUE, withSmooth = TRUE, df = NULL)
##
.computeInterpolators(myRDA0, myRDA,withSmoothFct = TRUE)
