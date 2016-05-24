# Testing code for the RCMIP5 (?) 'checkTimePeriod.R' script

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf
library(testthat)

# To run this code:
#   source("checkTimePeriod.R")
#   library(testthat)
#   test_file("tests/testthat/test_checkTimePeriod.R")

context("checkTimePeriod")

test_that("checkTimePeriod handles bad input", {
    path <- "../../sampledata"
    if(!file.exists(path)) skip("Path doesn't exist")
    
    d <- getFileInfo(path)
    expect_error(checkTimePeriod())                     # no input
    expect_error(checkTimePeriod(c(d, d)))               # multi-value input
    expect_error(checkTimePeriod(12))                   # non-data.frame input
    expect_error(checkTimePeriod(data.frame()))         # incorrect data.frame input
})

test_that("checkTimePeriod correctly finds missing files", {
    path <- ("testdata_missingfile")
    d <- checkTimePeriod(getFileInfo(path))
    expect_is(d,"data.frame")
    expect_equal(nrow(d), 2)     # should be two cases
    expect_equal(ncol(d), 10)
    expect_false(d$allHere[1])  # monthly data case is not complete
    expect_false(d$allHere[2])  # annual data case is not complete
    expect_true(all(d$files > 1))
})

test_that("checkTimePeriod correctly sees continuous files", {
    path <- "../../sampledata"
    if(!file.exists(path)) skip("Path doesn't exist")
    
    d <- checkTimePeriod(getFileInfo(path))
    expect_is(d, "data.frame")
    expect_more_than(nrow(d), 2)     # should be several cases
    expect_equal(ncol(d), 10)
    expect_true(all(d$allHere))
})

test_that("checkTimePeriod correctly parses dates", {
    path <- "../../sampledata"
    if(!file.exists(path)) skip("Path doesn't exist")
    
    d <- checkTimePeriod(getFileInfo(path))
    expect_is(d, "data.frame")
    expect_is(d$startDate, "numeric")
    expect_is(d$endDate, "numeric")
    expect_true(all(d$endDate >= d$startDate))
    expect_true(all(d$startDate >= 1850))
    expect_true(all(d$endDate <= 2300))
})

test_that('checkTimePeriod correctly flags sub-monthly ensembles', {
    path <- ('testdata_shortFreq')
    d <- checkTimePeriod(getFileInfo(path))
    expect_true(is.na(d$allHere))
})
