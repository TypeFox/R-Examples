# Testing code for the RCMIP5 'makeMonthlyStat.R' script

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf
library(testthat)

# To run this code: 
#   source("makeMonthlyStat.R")
#   source("RCMIP5.R") # for cmip5data
#   library(testthat)
#   test_file("tests/testthat/test_makeMonthlyStat.R")

context("makeMonthlyStat")

test_that("makeMonthlyStat handles bad input", {
    expect_error(makeMonthlyStat(1))                         # non-list d
    expect_error(makeMonthlyStat(cmpi5data()))               # wrong size list d
    expect_error(makeMonthlyStat(d,verbose=1))               # non-logical verbose
    expect_error(makeMonthlyStat(d,verbose=c(F, F)))          # multiple verbose values
    expect_error(makeMonthlyStat(d,FUN=1))                   # non-function FUN
    expect_error(makeMonthlyStat(d,FUN=c(mean, mean)))        # multiple FUN values
})

test_that("makeMonthlyStat handles monthly data", {
    years <- 1850:1851
    d <- cmip5data(years)
    res <- makeMonthlyStat(d, verbose=F)
    
    # Is 'res' correct type and size?
    expect_is(res,"cmip5data")
    
    # Did unchanging info get copied correctly?
    expect_equal(res$lon, d$lon)
    expect_equal(res$lat, d$lat)
    expect_equal(res$valUnit, d$valUnit)
    expect_equal(res$files, d$files)
    
    # numYears set and provenance updated?
    expect_is(res$numYears, "integer")
    expect_more_than(nrow(res$provenance), nrow(d$provenance))
    
    # Does time match what we expect?
    expect_equal(res$time, 1:12)
    
    # Is the answer value data frame correctly sized?
    expect_equal(nrow(res$val), length(d$lon) * length(d$lat) * 12) 
    expect_equal(length(res$time), 12)
    
    # Are the answer values numerically correct?    
    d$val$time <- floor((d$val$time %% 1) * 12 + 1)
    dummyans <- aggregate(value~lon+lat+time, data=d, FUN=mean)
    expect_equal(dummyans$value, res$val$value)
})

test_that("makeMonthlyStat handles annual data", {
    years <- 1850:1851
    d <- cmip5data(years, monthly=F)
    expect_error(makeMonthlyStat(d, verbose=F))
})

test_that("makeMonthlyStat handles 4-dimensional data", {
    years <- 1850:1851
    d <- cmip5data(years, Z=T)
    res <- makeMonthlyStat(d, verbose=F)
    
    # Do years match what we expect?
    expect_equal(res$time, 1:12)
    
    # Is the answer value data frame correctly sized?
    expect_equal(nrow(res$val), length(d$lon) * length(d$lat) * length(d$Z) * 12) 
    expect_equal(length(res$time), 12)
    
    # Are the answer values numerically correct?    
    d$val$time <- floor((d$val$time %% 1) * 12 + 1)
    dummyans <- aggregate(value~lon+lat+Z+time, data=d, FUN=mean)
    expect_equal(dummyans$value, res$val$value)
})

test_that("makeMonthlyStat handles custom function and dots", {
    years <- 1850:1851
    llsize <- 2
    d <- cmip5data(years, lonsize=llsize, latsize=llsize)
    
    # All 1850 data 1, all 1851 data 2
    d$val$value <- 1
    d$val$value[floor(d$val$time) == years[2]] <- 2    
    w <- c(3, 1)
    
    # Compute correct answer
    d$val$month <- floor((d$val$time %% 1) * 12 + 1)
    ans <- aggregate(value~lon+lat+month, data=d$val, FUN=weighted.mean, w=w)
    
    res1 <- makeMonthlyStat(d, verbose=F, FUN=weighted.mean, w)
    expect_is(res1, "cmip5data")
    
    myfunc <- function(x, w, ...) weighted.mean(x, w, ...)
    res2 <- makeMonthlyStat(d, verbose=F, FUN=myfunc, w)
    expect_is(res1, "cmip5data")
    
    # Are the answer values numerically correct?    
    expect_equal(res1$val$value, ans$value)    
    expect_equal(res2$val$value, ans$value)
})
