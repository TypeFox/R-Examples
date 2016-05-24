# Test sss integration
#
# Author: Andrie
#------------------------------------------------------------------------------


if(interactive()){
  library(testthat)
  sampleRoot <- "sss/tests/testthat/samples/sample-1"
} else {
  sampleRoot <- "samples/sample-1"
}
filenameSSS <- file.path(sampleRoot, "recruitment_test.sss")
filenameASC <- file.path(sampleRoot, "recruitment_test.csv")
file.exists(filenameSSS)

#------------------------------------------------------------------------------

context("sample 1")

test_that("parsing of .sss and .asc works", {
	test <- read.sss(filenameSSS, filenameASC)
	expect_is(test, "data.frame")
})

