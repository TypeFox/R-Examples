pkg <- "rstpm2"
require(RUnit, quietly = TRUE)
library(package = pkg, character.only = TRUE)
## path <- system.file("unitTests", package = pkg)
path <- "c:/usr/src/R/rstpm2/inst/unitTests"
testSuite <- defineTestSuite(name=paste(pkg, "unit testing"),
                             dirs=path)

tests <- runTestSuite(testSuite)
