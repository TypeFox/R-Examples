if(require(testthat)) {
  require(traitr)

  context("itemlist")
  i <- itemList(items=list(),
                item_factory=function(.) numericItem(0),
                name="test")
  ## get_test is alias to get_value
  expect_that(length(i$get_test()) == 0, is_true())
  
  ## add item
  i$append_item(i$item_factory())
  expect_that(length(i$get_value()) == 1, is_true())

  ## remove item
  i$append_item(i$item_factory())
  i$remove_item(2)
  expect_that(length(i$get_value()) == 1, is_true())  

  ## to_R
  i$to_R <- function(.) sapply(i$get_value(), function(j) j$to_R())
  expect_that(i$to_R()[[1]] == 0, is_true())



}
