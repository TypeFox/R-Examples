# library(testthat)
# 
# dsLinear <- read.table(file="./inst/extdata/BirthRatesOk.txt", header=TRUE, sep="\t", stringsAsFactors=F)
# dsLinear$Date <- as.Date(dsLinear$Date)
# changeMonth <- as.Date("1996-02-15") # as.Date(dateBombing + weeks(40))
# dsLinear$StageID <- ifelse(dsLinear$Date < changeMonth, 1L, 2L)
# 
# ###########
# context("Linear")
# ###########
# test_that("Smoke Test", {
# #   Wats::LinearPlot(dsPlot=dsLinear, xName="Date", yName="BirthRate", idName="StageID")
#   
# #   expect_equal(returned_object$data, expected=data.frame(), label="An empty data.frame should be returned.")
# #   expect_equal(returned_object$raw_csv, expected=raw(0))
# #   expect_true(is.null(returned_object$records_collapsed))
# #   expect_true(is.null(returned_object$fields_collapsed))
# #   expect_equal(returned_object$status_message, expected="Reading the REDCap data was not successful.  The error message was:\nError in textConnection(text) : invalid 'text' argument\n")
# #   expect_false(returned_object$success)
# })
