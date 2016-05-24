# Testing code for the RCMIP5 'saveNetCDF.R' script

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf
library(testthat)

# To run this code: 
#   library(testthat)
#   library(ncdf4)
#   source("saveNetCDF.R")
#   source("RCMIP5.R") # for cmip5data
#   test_file("tests/testthat/test_saveNetCDF.R")

context("saveNetCDF")

test_that("saveNetCDF handles bad input", {
    if(!require(ncdf4, quietly=T)) skip("ncdf4 not available")
    
    d <- cmip5data(1)
    expect_error(saveNetCDF(1))                       # non-cmip5data x
    expect_error(saveNetCDF(d, file=1))               # non-character file
    expect_error(saveNetCDF(d, file=c("x", "x")))     # multiple file values
    expect_error(saveNetCDF(d, path=1))               # non-character path
    expect_error(saveNetCDF(d, path=c("x", "x")))     # multiple path values
    expect_error(saveNetCDF(d, verbose=1))            # non-logical verbose
    expect_error(saveNetCDF(d, verbose=c(T, T)))      # multiple verbose values
})

test_that("saveNetCDF saves X-Y-T data correctly", {
    if(!require(ncdf4, quietly=T)) skip("ncdf4 not available")
    
    d <- cmip5data(1)
    dfile <- tempfile()
    if(file.exists(dfile)) expect_true(file.remove(dfile))
    saveNetCDF(d, file=basename(dfile), path=dirname(dfile), verbose=F)
    expect_true(file.exists(dfile))
    
    if(file.exists(dfile)) {
        nc <- nc_open(dfile)
        test <- ncvar_get(nc, "var")
        
        expect_equal(length(nc$var), 1) # variable count should match
        expect_equal(length(nc$dim), 3) # dimension count should match
        expect_equivalent(d$val$value, as.numeric(ncvar_get(nc, "var")))  # data should match
        expect_is(nc$dim$lon$units, "character")  # units should be written...
        expect_is(nc$dim$lat$units, "character")
        expect_null(nc$dim$Z$units)
        expect_is(nc$dim$time$units, "character")
        expect_is(nc$dim$time$units, "character")
        expect_is(nc$dim$time$calendar, "character")
    }
})

test_that("saveNetCDF saves X-Y-Z-T data correctly", {
    if(!require(ncdf4, quietly=T)) skip("ncdf4 not available")
    
    d <- cmip5data(1, Z=T)
    dfile <- tempfile()
    if(file.exists(dfile)) expect_true(file.remove(dfile))
    saveNetCDF(d, file=basename(dfile), path=dirname(dfile), verbose=F)
    expect_true(file.exists(dfile))
    
    if(file.exists(dfile)) {
        nc <- nc_open(dfile)
        test <- ncvar_get(nc, "var")

        expect_equal(length(nc$var), 1) # variable count should match
        expect_equal(length(nc$dim), 4) # dimension count should match
        expect_equivalent(d$val$value, as.numeric(ncvar_get(nc, "var")))  # data should match
        expect_is(nc$dim$lon$units, "character")  # units should be written...
        expect_is(nc$dim$lat$units, "character")
        expect_is(nc$dim$Z$units, "character")
        expect_is(nc$dim$time$units, "character")
        expect_is(nc$dim$time$units, "character")
        expect_is(nc$dim$time$calendar, "character")
    } 
})

test_that("saveNetCDF saves X-Y (area) data correctly", {
    if(!require(ncdf4, quietly=T)) skip("ncdf4 not available")
    
    d <- cmip5data(0, time=F)  # area data
    dfile <- tempfile()
    if(file.exists(dfile)) expect_true(file.remove(dfile))
    saveNetCDF(d, file=basename(dfile), path=dirname(dfile), verbose=F)
    expect_true(file.exists(dfile))
    
    if(file.exists(dfile)) {
        nc <- nc_open(dfile)
        test <- ncvar_get(nc)
        
        expect_equal(length(nc$var), 1) # variable count should match
        expect_equal(length(nc$dim), 2) # dimension count should match
        expect_equivalent(d$val$value, as.numeric(ncvar_get(nc)))  # data should match
        expect_is(nc$dim$lon$units, "character")  # units should be written...
        expect_is(nc$dim$lat$units, "character")
        expect_null(nc$dim$depth$units)
        expect_null(nc$dim$lev$units)
        expect_null(nc$dim$time$units)
        expect_null(nc$dim$time$calendar)
    }
})

test_that("saveNetCDF saves t (time) data correctly", {
    if(!require(ncdf4, quietly=T)) skip("ncdf4 not available")
    
    d <- cmip5data(1:10, lonlat=F)  # time data
    dfile <- tempfile()
    if(file.exists(dfile)) expect_true(file.remove(dfile))
    saveNetCDF(d, file=basename(dfile), path=dirname(dfile), verbose=F)
    expect_true(file.exists(dfile))
    
    if(file.exists(dfile)) {
        nc <- nc_open(dfile)
        test <- ncvar_get(nc)
        
        expect_equal(length(nc$var), 1) # variable count should match
        expect_equal(length(nc$dim), 1) # dimension count should match
        expect_equivalent(d$val$value, as.numeric(ncvar_get(nc)))  # data should match
        expect_null(nc$dim$lon)
        expect_null(nc$dim$lat)
        expect_null(nc$dim$Z)
        expect_null(nc$dim$depth$units)
        expect_null(nc$dim$lev$units)
        expect_is(nc$dim$time$units, "character")
        expect_is(nc$dim$time$calendar, "character")
    }
})
