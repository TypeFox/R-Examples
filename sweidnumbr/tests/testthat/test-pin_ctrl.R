
context("pin_ctrl")

pin_test <- c("196408233234", "196408233235")
pin_test_res <- c(TRUE, FALSE)

test_that(desc="control number",{
  expect_is(pin_ctrl(pin = pin_test), "logical")
  expect_equal(pin_ctrl(pin = pin_test), expected = pin_test_res)
})

test_that(desc="Handle NA, interim and coordn in pin_ctrl",{
  suppressWarnings(expect_equal(pin_ctrl(pin = c("195812793098", "195812793099", "19581279P092", "19581279P091")),
               c(TRUE, FALSE, NA, NA)))
  expect_true(is.na(pin_ctrl(c(NA,"198501169885"))[1]))
  expect_false(is.na(pin_ctrl(c(NA,"198501169885"))[2]))
})

test_that(desc="Expect force_logical",{
  num_to_check <- c("202100-6255","121212-1212","19121212-1212","121212+1212","1212121212",
                    1212121212,"19121212+1212", NA, Inf, TRUE, F, "foo", 123, 456L)
  suppressWarnings(expect_equal(pin_ctrl(num_to_check), 
               c(NA, TRUE, TRUE, TRUE, TRUE, TRUE, NA, NA, NA, NA, NA, NA, NA, NA)))
  expect_equal(pin_ctrl(num_to_check, force_logical=TRUE), 
               c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
})
