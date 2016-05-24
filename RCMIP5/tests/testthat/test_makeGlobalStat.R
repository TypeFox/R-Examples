# Testing code for the RCMIP5 'makeGlobalStat.R' script

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf
library(testthat)

# To run this code: 
#   source("makeGlobalStat.R")
#   source("RCMIP5.R") # for cmip5data
#   library(testthat)
#   test_file("tests/testthat/test_makeGlobalStat.R")

context("makeGlobalStat")

test_that("makeGlobalStat handles bad input", {
    expect_error(makeGlobalStat(1))                         # non-list d
    expect_error(makeGlobalStat(cmpi5data()))               # wrong size list d
    expect_error(makeGlobalStat(d,verbose=1))               # non-logical verbose
    expect_error(makeGlobalStat(d,verbose=c(T, T)))          # multiple verbose values
    expect_error(makeGlobalStat(d,FUN=1))                   # non-function FUN
    expect_error(makeGlobalStat(d,FUN=c(mean, mean)))        # multiple FUN values
})

test_that("makeGlobalStat handles monthly data", {
    years <- 1850:1851
    d <- cmip5data(years)
    res <- makeGlobalStat(d, verbose=F)
    
    # Is 'res' correct type and size?
    expect_is(res,"cmip5data")
    
    # Did unchanging info get copied correctly?
    expect_equal(res$valUnit, d$valUnit)
    expect_equal(res$files, d$files)
    
    # Lon/lat removed, numCells set, and provenance updated?
    expect_null(res$lon)
    expect_null(res$lat)
    expect_is(res$numCells, "integer")
    expect_more_than(nrow(res$provenance), nrow(d$provenance))
    
    # Does time match what we expect?
    expect_equal(res$time, d$time)
    
    # Is the answer value data frame correctly sized?
    expect_equal(nrow(res$val), nrow(d$val)/length(d$lon)/length(d$lat))
    
    # Are the answer values numerically correct?
    expect_equal(mean(res$val$value), mean(d$val$value))  # no weighting
})

test_that("makeGlobalStat weights correctly", {
    d <- cmip5data(1850, randomize=T, monthly=F)
    darea <- cmip5data(0, time=F, randomize=T)  # create an area file
    
    res <- makeGlobalStat(d, area=darea, verbose=F)
    
    # Are the answer values numerically correct?
    dummyans <- weighted.mean(d$val$value, w=darea$val$value)
    expect_equal(dummyans, res$val$value)
})

test_that("weighted.sum works correctly", {
    d <- cmip5data(1850, randomize=T, monthly=F)
    darea <- cmip5data(0, time=F, randomize=T)
    res <- makeGlobalStat(d, area=darea, verbose=F, FUN=weighted.sum)
    
    # Are the answer values numerically correct?
    dummyans <- weighted.sum(d$val$value, w=darea$val$value)
    expect_equal(dummyans, res$val$value)
    
    # Make sure the function itself is OK
    expect_equal(weighted.sum(1:4), 10)
    expect_equal(weighted.sum(1:4, 1:4), 30) # 4*4 + 3*3 + 2*2 + 1*1
})

test_that("makeGlobalStat handles 4-dimensional data", {
    years <- 1850:1851
    d <- cmip5data(years, Z=T)
    res <- makeGlobalStat(d, verbose=F)
    
    # Do years match what we expect?
    expect_equal(res$time, d$time)
    
    # Is the answer value array correctly sized?
    expect_equal(nrow(res$val), nrow(d$val)/length(d$lon)/length(d$lat))
})

test_that("makeGlobalStat handles custom function and dots", {
    years <- 1850:1851
    llsize <- 2
    d <- cmip5data(years, lonsize=llsize, latsize=llsize)
    darea <- cmip5data(0, time=F, lonsize=llsize, latsize=llsize)
    
    # All data 1 except for max lon/lat is 2
    d$val$value <- 1
    d$val$value[d$val$lon == max(d$lon) & d$val$lat == max(d$lat)] <- 2    
    darea$val$value <- c(rep(1, llsize*llsize-1), llsize*llsize-1)
    
    # Compute correct answer
    ans <- aggregate(value~time, data=d$val, FUN=weighted.mean, w=darea$val$value)
        
    res1 <- makeGlobalStat(d, darea, verbose=F, FUN=weighted.mean)
    expect_is(res1, "cmip5data")
    
    myfunc <- function(x, w, ...) weighted.mean(x, w, ...)
    res2 <- makeGlobalStat(d, darea, verbose=F, FUN=myfunc)
    expect_is(res1, "cmip5data")
    
    # Are the result values correct?    
    expect_equal(res1$val$value, ans$value)    
    expect_equal(res2$val$value, ans$value)
})

test_that("makeGlobalStat sorts before computing", {
    years <- 1850:1851
    llsize <- 2
    d <- cmip5data(years, lonsize=llsize, latsize=llsize, monthly=F)
    darea <- cmip5data(0, time=F, lonsize=llsize, latsize=llsize)
    
    # All data 1 except for max lon/lat is 2
    d$val$value <- 1
    d$val$value[d$val$lon == max(d$lon) & d$val$lat == max(d$lat)] <- 2    
    darea$val$value <- c(rep(1, llsize*llsize-1), llsize*llsize-1)
    
    # Compute correct answer
    ans <- aggregate(value~time, data=d$val, FUN=weighted.mean, w=darea$val$value)
    
    # Now we put `darea` out of order and call makeGlobalStat
    darea$val <- arrange(darea$val, desc(lon), desc(lat))
    res1 <- makeGlobalStat(d, darea, verbose=F)
    expect_is(res1, "cmip5data")

    # makeGlobalStat should be sorted darea correctly before calculating
    # Are the result values correct?    
    expect_equal(res1$val$value, ans$value)
    
    # Put data out of order and test again
    d$val <- arrange(d$val, desc(lon), desc(lat))
    res2 <- makeGlobalStat(d, darea, verbose=F)
    expect_equal(res2$val$value, ans$value)
})


