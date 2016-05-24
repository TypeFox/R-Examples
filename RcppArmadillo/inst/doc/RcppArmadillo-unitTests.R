### R code from vignette source 'RcppArmadillo-unitTests.Rnw'

###################################################
### code chunk number 1: RcppArmadillo-unitTests.Rnw:12-16
###################################################
require(RcppArmadillo)
prettyVersion <- packageDescription("RcppArmadillo")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")
library(RUnit)


###################################################
### code chunk number 2: unitTesting
###################################################
pkg <- "RcppArmadillo"
if (file.exists("unitTests-results")) unlink("unitTests-results", recursive = TRUE)
dir.create("unitTests-results")
pathRcppArmadilloTests <<- system.file("unitTests", package = pkg)
path <- system.file("unitTests", package=pkg)
testSuite <- defineTestSuite(name=paste(pkg, "unit testing"), dirs=path)
tests <- runTestSuite(testSuite)
err <- getErrors(tests)
if (err$nFail > 0) stop(sprintf("unit test problems: %d failures", err$nFail))
if (err$nErr > 0) stop( sprintf("unit test problems: %d errors", err$nErr))
printHTMLProtocol(tests, fileName=sprintf("unitTests-results/%s-unitTests.html", pkg))
printTextProtocol(tests, fileName=sprintf("unitTests-results/%s-unitTests.txt" , pkg))

#if (file.exists("/tmp")) {
#    invisible(sapply(c("txt", "html"), function(ext) {
#        fname <- sprintf("unitTests-results/%s-unitTests.%s", pkg, ext)
#        file.copy(fname, "/tmp", overwrite=TRUE)
#    }))
#}


###################################################
### code chunk number 3: importResults
###################################################
results <- "unitTests-results/RcppArmadillo-unitTests.txt"
if (file.exists(results)) {
    writeLines(readLines(results))
} else{
    writeLines( "unit test results not available" )
}


