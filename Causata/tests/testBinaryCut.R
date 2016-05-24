library(testthat)
library(Causata)

#source("utils.R")
equals <- testthat::equals

context("BinaryCut")

x1    <- c(1,1,1,1, 2,2,2,2,2, 3,3,3,3,3,3,  4, 4, 4)
x2    <- c(1,1,1,1, 2,2,2,2,2, 3,3,3,3,3,3, NA,NA,NA) # includes NAs
dv.a  <- c(1,1,1,0, 1,0,0,0,0, 1,1,1,0,0,0,  1, 1, 0) # last two levels have similar odds
dv.b  <- c(1,1,1,0, 1,0,0,0,0, 1,0,0,0,0,0,  1, 1, 0) # middle two levels have similar odds
breaks1  <- c(0,1,2,3,4)
breaks2  <- c(0,1,2,3)
woe1.a <- Woe(cut(x1, breaks1, include.lowest=TRUE), dv.a)
woe1.b <- Woe(cut(x1, breaks1, include.lowest=TRUE), dv.b)

bcl1.a        <- BinaryCut(x1, dv.a, bins=TRUE) # defaults, no missing
test_that("Five breaks returned with default values, no missing", expect_equal(bcl1.a$breaks, breaks1))

bcl2.a        <- BinaryCut(x2, dv.a, bins=TRUE) # defaults with missing
test_that("Four breaks returned with default values, includes missing", expect_equal(bcl2.a$breaks, breaks2))

f1.a.nbin2    <- BinaryCut(x1, dv.a, nbins=2) # only 2 bins
test_that("Setting nbins=2 returns 2 bins.", expect_equal(2, length(levels(f1.a.nbin2))))

f1.a.woemerge <- BinaryCut(x1, dv.a, woeDelta=1)
test_that("Last two levels with similar WOE are merged using woeDelta=1", expect_equal("(2,4]", levels(f1.a.woemerge)[3]))

f1.b.woemerge <- BinaryCut(x1, dv.b, woeDelta=1)
test_that("Middle two levels with similar WOE are merged using woeDelta=1", expect_equal("(1,3]", levels(f1.b.woemerge)[2]))

f1.a.nmerge   <- BinaryCut(x1, dv.a, minBin=2)
test_that("Smaller end levels are merged using minBin=2", expect_equal(c("[0,2]","(2,4]"), levels(f1.a.nmerge)))
