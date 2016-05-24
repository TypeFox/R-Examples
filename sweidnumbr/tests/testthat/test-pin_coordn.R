
context("pin_coordn")

pin_test <- c("196408233234", "196408830000")
pin_test_res <- c(FALSE, TRUE)

test_that(desc="control number",{
  expect_is(pin_coordn(pin = pin_test), "logical")
  expect_equal(pin_coordn(pin = pin_test), expected = pin_test_res)
})

test_that(desc="Handle NA",{
  expect_true(is.na(pin_coordn(as.pin(c(NA,"198501169885")))[1]))
  expect_false(is.na(pin_coordn(as.pin(c(NA,"198501169885")))[2]))
})
