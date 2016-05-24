# Testing code for the RCMIP5 'filterDimensions.R' script

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf
library(testthat)

# To run this code: 
#   library(testthat)
#   source("filterDimensions.R")
#   source("RCMIP5.R") # for cmip5data
#   test_file("tests/testthat/test_filterDimensions.R")

context("filterDimensions")

test_that("filterDimensions handles bad input", {
    expect_error(filterDimensions(1))                       # non-cmip5data x
    d <- cmip5data(1)
    expect_error(filterDimensions(d, lons='1'))              # non-numeric lons
    expect_error(filterDimensions(d, lats='1'))              # non-numeric lats
    expect_error(filterDimensions(d, Zs='1'))              # non-numeric Zs
    expect_error(filterDimensions(d, years='1'))              # non-numeric years
    expect_error(filterDimensions(d, months='1'))              # non-numeric months
    expect_error(filterDimensions(d, verbose=1))            # non-logical verbose
    expect_error(filterDimensions(d, verbose=c(T, T)))      # multiple verbose values
})

test_that("filterDimensions filters lon", {
    d <- cmip5data(1)
    d$lon <- NULL
    expect_warning(filterDimensions(d, lons=1))
    
    d <- cmip5data(1)
    lf <- d$lon[1:(length(d$lon)-1)] # the filter
    res <- filterDimensions(d, lons=lf)
    expect_equal(res$lon, lf)
    expect_equal(nrow(res$val), prod(length(res$lon), length(res$lat), length(res$time)))
    expect_more_than(nrow(res$provenance), nrow(d$provenance))     
})

test_that("filterDimensions filters lat", {
    d <- cmip5data(1)
    d$lat <- NULL
    expect_warning(filterDimensions(d, lats=1))
    
    d <- cmip5data(1)
    lf <- d$lat[-length(d$lat)] # the filter
    res <- filterDimensions(d, lats=lf)
    expect_equal(res$lat, lf)
    expect_equal(nrow(res$val), prod(length(res$lon), length(res$lat), length(res$time)))
    expect_more_than(nrow(res$provenance), nrow(d$provenance))     
})

test_that("filterDimensions filters Z", {
    expect_warning(filterDimensions(cmip5data(1), Zs=1))
    
    d <- cmip5data(1, Z=T)
    zf <- d$Z[-length(d$Z)] # the filter
    res <- filterDimensions(d, Zs=zf)
    expect_equal(res$Z, zf)
    expect_equal(nrow(res$val), prod(length(res$lon), length(res$lat), length(res$Z), length(res$time)))
    expect_more_than(nrow(res$provenance), nrow(d$provenance))     
})

test_that("filterDimensions filters time (years)", {
    
    # Annual data
    d <- cmip5data(1:5, monthly=F)
    d$time <- NULL
    expect_warning(filterDimensions(d, years=1))
    
    d <- cmip5data(1:5, monthly=F)
    yf <- d$time[-length(d$time)] # the filter
    res <- filterDimensions(d, years=yf)
    expect_equal(res$time, yf)                                      # only years in filter
    expect_equal(nrow(res$val), prod(length(res$lon), length(res$lat), length(res$time)))
    expect_more_than(nrow(res$provenance), nrow(d$provenance))  # provenance bigger
    
    # Monthly data
    d <- cmip5data(1:5)
    yf <- d$time[1:(length(d$time)/2)] # the filter
    res <- filterDimensions(d, years=yf)
    expect_equal(unique(floor(res$time)), unique(floor(yf)))        # only years in filter
    expect_equal(nrow(res$val), prod(length(res$lon), length(res$lat), length(res$time)))
    expect_more_than(nrow(res$provenance), nrow(d$provenance))  # and provenance bigger
})

test_that("filterDimensions filters time (months)", {
    
    # Annual data
    d <- cmip5data(1:5, monthly=F)
    expect_warning(filterDimensions(d, months=1))                   # monthly data required
    
    # Monthly data
    d <- cmip5data(1:5)
    mf <- 1:2
    res <- filterDimensions(d, months=mf)
    expect_equal(length(unique(round(res$time %% 1, 2))), length(mf))  # only months in filter
    expect_equal(nrow(res$val), prod(length(res$lon), length(res$lat), length(res$time)))
    expect_more_than(nrow(res$provenance), nrow(d$provenance))  # and provenance bigger
    
    expect_error(filterDimensions(d, months=13))                    # illegal months value
})

test_that("filterDimensions handles multiple operations", {
    d <- cmip5data(1:5, Z=T)
    res <- filterDimensions(d, lons=d$lon[1], lats=d$lat[1], Zs=d$Z[1], years=d$time[1])
    expect_equal(nrow(res$val), prod(length(res$lon), length(res$lat), length(res$time)))    
})

test_that("filterDimensions handles nonstandard structures", {
    d <- cmip5data(0, time=F)  # area-only data
    res <- filterDimensions(d, lons=d$lon[1], lats=d$lat[1])
    expect_is(res, "cmip5data")
    
    d <- cmip5data(0, Z=T, time=F)  # area and Z-only data
    res <- filterDimensions(d, lons=d$lon[1], lats=d$lat[1], Zs=d$Z[1])
    expect_is(res, "cmip5data")
    
    d <- cmip5data(1:5, lonlat=F)  # time-only data
    res <- filterDimensions(d, years=d$time[1])
    expect_is(res, "cmip5data")    
})

