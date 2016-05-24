# Testing code for the RCMIP5 scripts in 'RCMIP5.R'

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf
library(testthat)

# To run this code:
#   source("RCMIP5.R")
#   library(testthat)
#   test_file("tests/testthat/test_RCMIP5.R")

context("RCMIP5")

test_that("cmip5data print method works", {
    d <- cmip5data(2000:2005)
    expect_output(print(d), "CMIP5")
    # not sure what to test here, except that no error
})

test_that("cmip5data summary method detects summaries", {
    d <- cmip5data(2000:2005, Z=T)
    expect_output(print(summary(d)), "CMIP5")
    
    # Summary should let user know data have been run through stat fn
    da <- makeAnnualStat(d)
    expect_output(print(summary(da)), "annual summary")
    dm <- makeMonthlyStat(d)
    expect_output(print(summary(dm)), "monthly summary")
    dz <- makeZStat(d)
    expect_output(print(summary(dz)), "Z summary")
    dg <- makeGlobalStat(d)
    expect_output(print(summary(dg)), "spatial summary")
    
    # Multiple stat functions should be detected
    dag <- makeGlobalStat(da)
    expect_output(print(summary(dag)), "annual summary")
    expect_output(print(summary(dag)), "spatial summary")
    daz <- makeZStat(da)
    expect_output(print(summary(daz)), "annual summary")
    expect_output(print(summary(daz)), "Z summary")
    dmg <- makeGlobalStat(dm)
    expect_output(print(summary(dmg)), "monthly summary")
    expect_output(print(summary(dmg)), "spatial summary")
    dmz <- makeZStat(dm)
    expect_output(print(summary(dmz)), "monthly summary")
    expect_output(print(summary(dmz)), "Z summary")
    
    # All filter functions should be detected
    df <- filterDimensions(d, years=2000:2002)
    expect_output(print(summary(df)), "filtered")
    df <- filterDimensions(d, months=1:6)
    expect_output(print(summary(df)), "filtered")
    df <- filterDimensions(d, lons=d$lon)
    expect_output(print(summary(df)), "filtered")
    df <- filterDimensions(d, lats=d$lat)
    expect_output(print(summary(df)), "filtered")
    df <- filterDimensions(d, Zs=d$Z)
    expect_output(print(summary(df)), "filtered")
})

test_that("as.data.frame works", {
    df <- as.data.frame(cmip5data(2000:2002, Z=T))
    expect_is(df, "data.frame")
    expect_equal(names(df), c("lon", "lat", "Z", "time", "value"))
})

test_that("as.array works", {
    arr <- as.array(cmip5data(2000:2002, Z=T))
    expect_is(arr, "array")
    expect_equal(dim(arr), c(10, 10, 5, 36))

    arr <- as.array(cmip5data(2000:2002))
    expect_is(arr, "array")
    expect_equal(dim(arr), c(10, 10, 36))
    
    arr <- as.array(cmip5data(2000:2002), drop=FALSE)
    expect_is(arr, "array")
    expect_equal(dim(arr), c(10, 10, 1, 36))
})
