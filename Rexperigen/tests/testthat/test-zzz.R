context("zzz")
library(Rexperigen)
## TODO: Rename context
## TODO: Add more tests

test_that("initialization is okay", {
  expect_equal(getOption("Rexperigen.experimenter"), "")
  expect_equal(getOption("Rexperigen.password"), "")
  expect_equal(getOption("Rexperigen.server"), "db.phonologist.org")
  expect_equal(getOption("Rexperigen.server.version"), "1.0.0")
})
