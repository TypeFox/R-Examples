### R code from vignette source 'RProtoBuf-unitTests.Rnw'

###################################################
### code chunk number 1: RProtoBuf-unitTests.Rnw:9-13
###################################################
require(RProtoBuf)
prettyVersion <- packageDescription("RProtoBuf")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")
library(RUnit)


###################################################
### code chunk number 2: unitTesting
###################################################
pkg <- "RProtoBuf"

if (file.exists("unitTests-results")) unlink("unitTests-results", recursive = TRUE)
dir.create("unitTests-results")
pathRcppTests <<- system.file("unitTests", package = pkg)
path <- system.file("unitTests", package=pkg)
testSuite <- defineTestSuite(name=paste(pkg, "unit testing"), dirs=path)
tests <- runTestSuite(testSuite)
err <- getErrors(tests)
if (err$nFail > 0) cat(sprintf("unit test problems: %d failures", err$nFail))
if (err$nErr > 0) cat(sprintf("unit test problems: %d errors", err$nErr))
printHTMLProtocol(tests, fileName=sprintf("unitTests-results/%s-unitTests.html", pkg))
printTextProtocol(tests, fileName=sprintf("unitTests-results/%s-unitTests.txt" , pkg))

if (file.exists("/tmp")) {
    invisible(sapply(c("txt", "html"), function(ext) {
        fname <- sprintf("unitTests-results/%s-unitTests.%s", pkg, ext)
        file.copy(fname, "/tmp", overwrite=TRUE)
    }))
}


###################################################
### code chunk number 3: importResults
###################################################
results <- sprintf("unitTests-results/%s-unitTests.txt", pkg)
if (file.exists(results)) {
    writeLines(readLines(results))
} else{
    writeLines("Unit test results not available")
}


