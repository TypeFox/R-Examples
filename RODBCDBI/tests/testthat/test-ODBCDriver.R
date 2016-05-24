context("ODBCDriver")

test_that("ODBCDriver had nothing about DB Info", {
  driver <- RODBCDBI::ODBC()
  expect_null(dbGetInfo(driver))
})
