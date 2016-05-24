context("yeastgfp")

test_that("yeastgfp() makes specific errors and messages", {
  expect_error(yeastgfp("aaa"), "aaa is not one of the subcellular locations")
  expect_message(yeastgfp("bud"), "no exclusive localization found for bud")
})

