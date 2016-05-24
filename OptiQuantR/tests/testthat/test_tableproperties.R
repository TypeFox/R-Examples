library(data.table)
library(testthat)
library(OptiQuantR)
context("Table properties")

# testdata <- system.file("extdata", "testdata.csv", package = "OptiQuantR")

names <- c("timeStampO",
            "session_id",
            "session_started",
            "session_finished",
            "session_total_counted",
            "gate_id",
            "person_kind",
            "id",
            "year",
            "month",
            "day",
            "weekday",
            "hour",
            "week",
            "min_start_off",
            "session_length",
            "timeStamp")

test_that("output has 17 columns", {
  expect_equal(ncol(readoqcsv(system.file("extdata", "testdata.csv", package = "OptiQuantR"))), 17)
})

test_that("output has the right column names", {
  expect_equal(colnames(readoqcsv(system.file("extdata", "testdata.csv", package = "OptiQuantR"))),
               names)
})

test_that("output is of class data.frame & data.table", {
  expect_true(is.data.frame(readoqcsv(system.file("extdata", "testdata.csv", package = "OptiQuantR"))))
  expect_true(is.data.table(readoqcsv(system.file("extdata", "testdata.csv", package = "OptiQuantR"))))
})

