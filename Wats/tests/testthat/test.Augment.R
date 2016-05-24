library(testthat)
#filePathOutcomes <- file.path(base::path.package("Wats"), "extdata", "BirthRatesOk.txt")
# filePathOutcomes <- file.path(devtools::inst(name="Wats"), "extdata", "BirthRatesOk.txt") #This approach accounts for working on developmental box.
# filePathOutcomes <- file.path(devtools::inst(name="Wats"), "extdata", "BirthRatesOk.txt") #This approach accounts for working on developmental box.


###########
context("Augment")
###########
test_that("AugmentYearDataWithMonthResolution", {
  dsBirthRate <- CountyMonthBirthRate2005Version[CountyMonthBirthRate2005Version$CountyName=="oklahoma", ]
#   dsBirthRate$Date <- as.Date(dsBirthRate$Date)
#   changeMonth <- as.Date("1996-02-15") # as.Date(dateBombing + weeks(40))
#   dsBirthRate$StageID <- ifelse(dsBirthRate$Date < changeMonth, 1L, 2L)
  
  expectedColumnNames <- c(colnames(dsBirthRate), "CycleTally", "ProportionThroughCycle", "ProportionID"
                           , "StartingPointInCycle", "TerminalPointInCycle", "StageProgress")
  actual <- Wats::AugmentYearDataWithMonthResolution(ds=dsBirthRate, dateName="Date")
  
  expect_equal(mean(actual$ProportionThroughCycle), expected=0.5, label="The mean of the proportion should be 0.5.")
  expect_equal(actual$ProportionID, expected=rep(1:12, times=10), label="The `ProportionID` should be correct.")
  expect_equal(actual$CycleTally, expected=rep(0:9, each=12), label="There should be 120 Cycle Indicies.")
  expect_equal(actual$StartingPointInCycle, expected=rep(c(T,F,F,F,F,F,F,F,F,F,F,F), time=10), label="The `StartingPointInCycle` should be correct.")
  expect_equal(actual$TerminalPointInCycle, expected=rep(c(F,F,F,F,F,F,F,F,F,F,F,T), time=10), label="The `TerminalPointInCycle` should be correct.")
  expect_equal(colnames(actual), expected=expectedColumnNames, label="The correct columns should be added.")
})

test_that("AugmentYearDataWithSecondResolution", {
  dsBirthRate <- CountyMonthBirthRate2005Version
  dsBirthRate <- dsBirthRate[dsBirthRate$CountyName=="oklahoma", ]
  dsBirthRate$Date <- as.POSIXct(dsBirthRate$Date, tz="GMT")
  
  expectedColumnNames <- c(colnames(dsBirthRate), "CycleTally", "ProportionThroughCycle", "ProportionID"
                           , "StartingPointInCycle", "TerminalPointInCycle", "StageProgress")
  actual <- Wats::AugmentYearDataWithSecondResolution(ds=dsBirthRate, dateName="Date")
  
  expect_equal(mean(actual$ProportionThroughCycle), expected=0.4933366, tolerance=1e-7, label="The mean of the proportion should be a little less than 0.5 (ie, ~.4933366) because the calender's first months are shorter than its last.")
  expect_equal(actual$ProportionID, expected=rep(1:12, times=10), label="The `ProportionID` should be correct.")
  expect_equal(actual$CycleTally, expected=rep(0:9, each=12), label="There should be 120 Cycle Indicies.")
  expect_equal(actual$StartingPointInCycle, expected=rep(c(T,F,F,F,F,F,F,F,F,F,F,F), time=10), label="The `StartingPointInCycle` should be correct.")
  expect_equal(actual$TerminalPointInCycle, expected=rep(c(F,F,F,F,F,F,F,F,F,F,F,T), time=10), label="The `TerminalPointInCycle` should be correct.")
  expect_equal(colnames(actual), expected=expectedColumnNames, label="The correct columns should be added.")
})
