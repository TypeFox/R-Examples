library("testthat")
library("permute")

context("Testing get methods")

test_that("default methods for get functions", {
    v <- 1:10
    expect_error(getBlocks(v), regexp = "No default method")
    expect_error(getPlots(v), regexp = "No default method")
    expect_error(getWithin(v), regexp = "No default method")
    expect_error(getStrata(v), regexp = "No default method")
    expect_error(getType(v), regexp = "No default method")
    expect_error(getMirror(v), regexp = "No default method")
    expect_error(getConstant(v), regexp = "No default method")
    expect_error(getNperm(v), regexp = "No default method")
    expect_error(getMaxperm(v), regexp = "No default method")
    expect_error(getMinperm(v), regexp = "No default method")
    expect_error(getMake(v), regexp = "No default method")
    expect_error(getObserved(v), regexp = "No default method")
    expect_error(getAllperms(v), regexp = "No default method")
    expect_error(getComplete(v), regexp = "No default method")
})
