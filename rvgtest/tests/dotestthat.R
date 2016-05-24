
## Location of unitTest files
runit.dir <- system.file(package="rvgtest", "unittests")

## Name of R file with test suite
testsuite.file <- "runalltests.R"

## path to R file with test suite
testsuite <- file.path(runit.dir, testsuite.file)

## Indicate that we run 'R CMD check' 
.rvgt.RCMDCHECK <- TRUE

## Load and test suite
cat("Running testsuite",testsuite,"\n") 
source(testsuite)

## End
                     
