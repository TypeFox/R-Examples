if(require(testthat)) {
  require(traitr)

  context("base")

  i <- BaseTrait$proto()

  ## has_slot
  expect_that(i$has_slot("a"), is_false())
  i$set_slot("a","b")
  expect_that(i$has_slot("a"), is_true())
  expect_that(i$get_slot("a") == "b", is_true())
  
  ## assign_if_null
  i$assign_if_null("test","value")
  expect_that(i$get_slot("test") == "value", is_true())
  i$assign_if_null("test","new value")
  expect_that(i$get_slot("test") == "value", is_true())

  ## do_call
  expect_that(is.null(i$do_call("test")), is_true())
  i$test <- function(.) 1
  expect_that(i$do_call("test") == 1, is_true())

  ## next_method
  a <- i$proto()
  a$test <- function(.) 2
  b <- a$proto()
  expect_that(b$next_method("test")(b) == 1, is_true())

  ## local slot
  a$a <- "b"
  expect_that(b$has_local_slot("a"), is_false())
  expect_that(a$has_local_slot("a"), is_true())
}
