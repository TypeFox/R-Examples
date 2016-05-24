# Testing code for the RCMIP5 scripts in 'RCMIP5.R'

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf
library(testthat)

# To run this code: 
#   source("RCMIP5.R")
#   library(testthat)
#   test_file("tests/testthat/test_internalHelpers.R")

context("cmip5data")

test_that("cmip5data handles bad input", {
    expect_error(cmip5data("hi"))   
    expect_error(cmip5data(1, monthly=123))   
    expect_error(cmip5data(1, Z=123))
    expect_error(cmip5data(1, lev=123))
    expect_error(cmip5data(1, randomize="hi"))   
})

test_that("cmip5data generates annual and monthly data", {
    d <- cmip5data(1)
    expect_is(d, "cmip5data")
    expect_equal(ncol(d$val), 5)
    expect_equal(all(d$val$Z), NA)
    expect_equal(length(unique(d$val$lon)), length(d$lon))
    expect_equal(length(unique(d$val$lat)), length(d$lat))
    expect_equal(length(unique(d$val$time)), length(d$time))    
    expect_equal(length(d$time), 12)
    expect_equal(d$debug$timeFreqStr, "mon")
    
    d <- cmip5data(1, monthly=F)
    expect_equal(ncol(d$val), 5)
    expect_equal(all(d$val$Z), NA)
    expect_equal(length(unique(d$val$lon)), length(d$lon))
    expect_equal(length(unique(d$val$lat)), length(d$lat))
    expect_equal(length(unique(d$val$time)), length(d$time))
    expect_equal(length(d$time), 1)
    expect_equal(d$debug$timeFreqStr, "yr")
})

test_that("cmip5data fills in ancillary data", {
    d <- cmip5data(1)
    expect_is(d$model, "character")
    expect_is(d$variable, "character")
    expect_is(d$experiment, "character")
    expect_is(d$valUnit, "character")
    expect_is(d$debug, "list")
    expect_is(d$debug$timeFreqStr, "character")
    expect_is(d$debug$calendarStr, "character")
    expect_is(d$debug$timeUnit, "character")
})

test_that("cmip5data obeys randomize", {
    expect_true(sum(cmip5data(1, randomize=T)$val$value) != sum(cmip5data(1, randomize=T)$val$value))
})

test_that("cmip5data creates area-only data", {
    lonsize <- 10
    latsize <- 10
    d <- cmip5data(1, lonsize=lonsize, latsize=latsize, time=F, verbose=F)
    expect_equal(ncol(d$val), 5)
    expect_equal(all(d$val$Z), NA)
    expect_equal(all(d$val$time), NA)
    expect_equal(length(unique(d$val$lon)), length(d$lon))
    expect_equal(length(unique(d$val$lat)), length(d$lat))
    expect_is(d$lon, "numeric")
    expect_is(d$lat, "numeric")
    expect_null(d$Z)
    expect_null(d$time)
})

test_that("cmip5data creates area and Z data", {
    lonsize <- 10
    latsize <- 10
    Zsize <- 5
    d <- cmip5data(1, lonsize=lonsize, latsize=latsize, Z=T, Zsize=Zsize, time=F, verbose=F)

    expect_equal(ncol(d$val), 5)
    expect_equal(all(d$val$time), NA)    
    expect_equal(length(unique(d$val$lon)), length(d$lon))
    expect_equal(length(unique(d$val$lat)), length(d$lat))
    expect_equal(length(unique(d$val$Z)), length(d$Z))
    expect_is(d$lon, "numeric")
    expect_is(d$lat, "numeric")
    expect_is(d$Z, "integer")
    expect_null(d$time) 
})

test_that("cmip5data creates time-only data", {
    d <- cmip5data(1, lonlat=F, Z=F, verbose=F)
    expect_equal(ncol(d$val), 5)
    expect_equal(all(d$val$lon), NA)
    expect_equal(all(d$val$lat), NA)
    expect_equal(all(d$val$Z), NA)
    expect_equal(length(unique(d$val$time)), length(d$time))
    expect_null(d$lon)
    expect_null(d$lat)
    expect_null(d$Z)
})
