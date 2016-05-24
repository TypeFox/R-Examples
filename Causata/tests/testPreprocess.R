library(testthat)
library(Causata)

context("PreProcess")

test_that("Converting R names to Causata names", {
  
  actual <- RToCausataNames(c("total.spend__All.Past", "page.view.count_Last.30.Days_1"))
  expect_equal(length(actual), 2)
  expect_equal(actual[["total.spend__All.Past"]], "total-spend$All Past")
  expect_equal(actual[["page.view.count_Last.30.Days_1"]], "page.view.count_Last.30.Days_1")
  
})

test_that("Converting Causata names to R names", {
  
  actual <- CausataToRNames(c("total-spend$All Past", "page-view-count$Last 30 Days"))
  expect_equal(length(actual), 2)
  expect_equal(actual[["total-spend$All Past"]], "total.spend__All.Past")
  expect_equal(actual[["page-view-count$Last 30 Days"]], "page.view.count__Last.30.Days")
  
})

test_that("ReplaceOutliers.CausataData works with strings", {

  df <- data.frame(variable1__All.Past=c(1,2,3,4,1000),
                   variable2__All.Past=c(5,6,7,8,1000) )
  causataData <- CausataData(df, rep(0,nrow(df)))
  causataData <- ReplaceOutliers(causataData, "variable1__All.Past", upperLimit=10)
  causataData <- ReplaceOutliers(causataData,  "variable2__All.Past" , upperLimit=11)
  expect_equal(causataData$variableList$variable1__All.Past$outlierUpperLimit, 10)
  expect_equal(causataData$variableList$variable2__All.Past$outlierUpperLimit, 11)

})