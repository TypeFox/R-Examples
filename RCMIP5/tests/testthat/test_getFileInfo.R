# Testing code for the RCMIP5 (?) 'getFileInfo.R' script

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf
library(testthat)

# To run this code:
#   source("getFileInfo.R")
#   library(testthat)
#   test_file("tests/testthat/test_getFileInfo.R")

context("getFileInfo")

test_that("getFileInfo handles bad input", {
    expect_error(getFileInfo("does_not_exist"),"file.exists\\(path) is not TRUE")
                                                            # path does not exist
    expect_error(getFileInfo(path=c("1","2")))              # multi-value path
    expect_error(getFileInfo(recursive=c(T,F)))             # multi-value recursive
    expect_error(getFileInfo(12))                           # non-character path
    expect_error(getFileInfo(recursive=1))                  # non-logical recursive
})

test_that("getFileInfo handles no input", {                 # no NetCDF files found
    w <- getOption('warn')
    options(warn=-1)
    path <- ("testdata_none")
    expect_warning(getFileInfo(path),"No NetCDF files found")
    expect_is(getFileInfo(path),"NULL")
    options(warn=w)
})

test_that("getFileInfo handles non-CMIP5 NetCDFs", {        # improper NetCDF filenames
    w <- getOption('warn')
    options(warn=-1)
    path <- ("testdata_badfilename")
    expect_warning(getFileInfo(path),"Unexpected")
    expect_is(getFileInfo(path),"NULL")
    options(warn=w)
})

test_that("getFileInfo handles annual NetCDFs", {
    path <- "../../sampledata/annual"
    if(!file.exists(path)) skip("Path doesn't exist")
    
    d <- getFileInfo(path)
    d <- d[complete.cases(d),]
    expect_is(d,"data.frame")
    expect_more_than(nrow(d),1)     # should be several
    expect_equal(ncol(d),9)
    expect_equal(d[1, "path"],path)
    expect_equal(d[1, "time"],"2171-2172")
})

test_that("getFileInfo handles monthly NetCDFs", {
    path <- "../../sampledata/monthly"
    if(!file.exists(path)) skip("Path doesn't exist")
           
    d <- getFileInfo(path)
    d <- d[complete.cases(d),]
    expect_is(d,"data.frame")
    expect_more_than(nrow(d),3)     # should be several
    expect_equal(ncol(d),9)
    expect_equal(d[1,"path"],path)
})

test_that("getFileInfo handles fixed NetCDFs", {
    path <- "../../sampledata/fx"
    if(!file.exists(path)) skip("Path doesn't exist")
         
    d <- getFileInfo(path)
    expect_is(d,"data.frame")
    expect_equal(ncol(d),9)
    expect_equal(d[1,"path"],path)
})

test_that("getFileInfo handles fixed and temporal NetCDFs together", {
    path <- "../../sampledata"
    if(!file.exists(path)) skip("Path doesn't exist")
        
    d <- getFileInfo(path)
    expect_is(d,"data.frame")
    expect_equal(ncol(d),9)
})

