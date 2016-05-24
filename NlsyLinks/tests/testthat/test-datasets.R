###########
context("Pair Datasets")
###########
test_that("Links79Pair", {
  expectedColumnCount <- 5
  actualColumnCount <- ncol(Links79Pair)
  expect_equal(actualColumnCount, expectedColumnCount, info="The number of columns should be correct.")
  
  expectedRowCount <- 42773 #42714 #11075
  actualRowCount <- nrow(Links79Pair)
  expect_equal(actualRowCount, expectedRowCount, info="The number of rows should be correct.")
})

test_that("Links79PairExpanded", {
  expectedColumnCount <- 26
  actualColumnCount <- ncol(Links79PairExpanded)
  expect_equal(actualColumnCount, expectedColumnCount, info="The number of columns should be correct.")
  
  expectedRowCount <- 42773 #42714 #11075
  actualRowCount <- nrow(Links79PairExpanded)
  expect_equal(actualRowCount, expectedRowCount, info="The number of rows should be correct.")
})
