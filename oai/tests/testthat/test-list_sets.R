context("list_sets")

test_that("list_sets", {
  skip_on_cran()

  aa <- list_sets()

  expect_is(aa, "data.frame")
  expect_is(aa, "oai_df")
  expect_is(aa$setSpec, "character")
  expect_is(aa$setName, "character")
})

# test_that("list_sets - resumption token", {
#   skip_on_cran()
#
#   aa <- list_sets(token = "1440473075578,null,null,50,null,null")
#   bb <- list_sets(token = "1440473094598,null,null,100,null,null")
#
#   expect_is(aa, "oai_df")
#   expect_is(bb, "oai_df")
#
#   expect_less_than(NROW(bb), NROW(aa))
# })

test_that("list_sets - curl options", {
  skip_on_cran()

  library("httr")

  expect_error(list_sets(config = timeout(0.001)), "Timeout was reached")
})

test_that("list_sets fails well", {
  skip_on_cran()

  expect_error(list_sets(token = 454),
               "The value of the resumptionToken argument is invalid or expired")
  expect_error(list_sets("stuff"),
               "One or more of your URLs")
})
