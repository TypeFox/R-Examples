
context("pin_sex")

pin_test <- c("198112189876", "196408233234", "198112180000", "187008233224", "000000")
pin_test_res <- as.factor(c("Male", "Male", "Female", "Female", NA))

test_that(desc="control number",{
  suppressWarnings(expect_is(pin_sex(pin = pin_test), "factor"))
  suppressWarnings(expect_equal(pin_sex(pin = pin_test), expected = pin_test_res))
})

test_that(desc="Handle NA, interim and coordn",{
  expect_true(is.na(pin_sex(c(NA,"198501169885"))[1]))
  expect_false(is.na(pin_sex(c(NA,"198501169885"))[2]))
  suppressWarnings(expect_equal(pin_sex(pin = c("19640823A234", "198112780000")), expected = as.factor(c(NA, "Female"))))
})
