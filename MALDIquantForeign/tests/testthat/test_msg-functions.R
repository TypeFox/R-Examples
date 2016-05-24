context("msg")


test_that("msg-functions", {
  expect_message(MALDIquantForeign:::.msg(TRUE, "foobar"), "foobar")
  expect_message(MALDIquantForeign:::.msg(TRUE, "foo", "bar"), "foobar")
  expect_that(MALDIquantForeign:::.msg(FALSE, "foobar"), not(shows_message()))
})
