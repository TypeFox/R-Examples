# Testing code for the RCMIP5 'mergeExperiments.R' script

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf
library(testthat)

# To run this code: 
#   source("mergeExperiments.R")
#   source("RCMIP5.R") # for cmip5data
#   library(testthat)
#   test_file("tests/testthat/test_mergeExperiments.R")

context("mergeExperiments")

test_that("mergeExperiments.R handles bad input", {    
    x <- cmip5data()
    expect_error(mergeExperiments())                             # non-cmip5data
    expect_error(mergeExperiments(1))                            # non-cmip5data
    expect_error(mergeExperiments(1, x))                         # non-cmip5data
    expect_error(mergeExperiments(x, 1))                         # non-cmip5data
    expect_error(mergeExperiments(x, x, verbose=1))              # non-logical verbose
    expect_error(mergeExperiments(x, x, verbose=c(T, T)))        # multiple verbose values
})

test_that("mergeExperiments identifies ancillary problems", {
    y <- cmip5data(2)
    
    x <- cmip5data(1)
    x$domain <- paste0(y$domain, "x")
    expect_error(mergeExperiments(x, y, verbose=F))
    
    x <- cmip5data(1)
    x$variable <- paste0(y$domain, "x")
    expect_error(mergeExperiments(x, y, verbose=F))
    
    x <- cmip5data(1)
    x$model <- paste0(y$domain, "x")
    expect_error(mergeExperiments(x, y, verbose=F))
    
    x <- cmip5data(1)
    x$valUnit <- paste0(y$valUnit, "x")
    expect_error(mergeExperiments(x, y, verbose=F))
    
    x <- cmip5data(1)
    x$lat <- c(y$lat, 0)
    expect_error(mergeExperiments(x, y, verbose=F))
    
    x <- cmip5data(1)
    x$lon <- c(y$lon, 0)
    expect_error(mergeExperiments(x, y, verbose=F))
    
    x <- cmip5data(1, Z=T)
    expect_error(mergeExperiments(x, y, verbose=F))
    
    x <- cmip5data(1)
    x$ensembles <- paste0(y$ensembles, "x")
    expect_warning(mergeExperiments(x, y, verbose=F))
})

test_that("mergeExperiments identifies time problems", {
    x <- cmip5data(1)
    y <- cmip5data(2)    
    x$debug$timeFreqStr <- paste0(y$debug$timeFreqStr, "x")
    expect_error(mergeExperiments(x, y, verbose=F))
 
    x <- cmip5data(1)
    expect_error(mergeExperiments(x, x, verbose=F)) # identical times
    x$time[length(x$time)] <- mean(y$time[1:2]) 
    expect_error(mergeExperiments(x, x, verbose=F)) # identical times
    
    x <- cmip5data(4)
    expect_warning(mergeExperiments(x, y, verbose=F)) # time gap
})

test_that("mergeExperiments merges monthly data", {
    x <- cmip5data(1:5)
    y <- cmip5data(6:10)    
    res <- mergeExperiments(x, y, verbose=F)
    
    expect_equal(nrow(x$val) + nrow(y$val), nrow(res$val))
    expect_equal(length(x$time) + length(y$time), length(res$time))
    expect_equal(c(x$files, y$files), res$files)
    expect_true(grepl(x$experiment, res$experiment)) # each experiment name should appear
    expect_true(grepl(y$experiment, res$experiment)) # each experiment name should appear
    
    # Order of x and y should make no difference
    xy <- mergeExperiments(x, y, verbose=F)
    yx <- mergeExperiments(x, y, verbose=F)
    expect_equal(xy$val, yx$val)
    expect_equal(xy$time, yx$time)
})

test_that("mergeExperiments merges annual data", {
    x <- cmip5data(1:5, monthly=F)
    y <- cmip5data(6:10, monthly=F)    
    res <- mergeExperiments(x, y, verbose=F)
    
    expect_equal(nrow(x$val) + nrow(y$val), nrow(res$val))
    expect_equal(length(x$time) + length(y$time), length(res$time))
})

test_that("mergeExperiments merges 4-dimensional data", {
    x <- cmip5data(1:5, Z=T)
    y <- cmip5data(6:10, Z=T)    
    res <- mergeExperiments(x, y, verbose=F)
    
    expect_equal(nrow(x$val) + nrow(y$val), nrow(res$val))
    expect_equal(length(x$time) + length(y$time), length(res$time))
})
