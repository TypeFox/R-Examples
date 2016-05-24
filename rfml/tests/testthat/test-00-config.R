context("Config")

test_that("MarkLogic database can be initiated for rfml", {
  skip_on_cran()
  ml.init.database(port = "8088")
})

test_that("MarkLogic database can be cleand up", {
  skip_on_cran()
  ml.clear.database(port = "8088")
  expect_error(ml.connect(port="8088"), "The database on http://localhost:8088 is not set up to work with rfml. Use ml.init.database for setting up the database.")
  ml.init.database(port = "8088")
})
