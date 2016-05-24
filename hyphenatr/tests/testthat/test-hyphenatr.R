context("basic functionality")
test_that("we can do something", {

  expect_that(hyphenate("testing"), equals("test=ing"))

  switch_dict("hu_HU")

  expect_that(curr_dict(), equals("hu_HU"))

  expect_that(length(list_dicts()), equals(39))

  expect_that(hyphendict_is_loaded(), equals(TRUE))

  switch_dict("de_DE")

  expect_that(hyphenate("tägelîch"), equals("tä=gelîch"))

})
