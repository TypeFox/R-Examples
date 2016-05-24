
context("pin_to_date")

test_that(desc="pin_to_date",{
  expect_equal(pin_to_date(pin = c("196408233234", "186408833224")), expected = lubridate::ymd(c("1964-08-23","1864-08-23")))
  check_class <- pin_to_date(pin = c("196408233234", "186408833224"))
  expect_true(inherits(check_class, "POSIXct") | inherits(check_class, "Date"))
})

test_that(desc="Handle NA and interimn in pin_to_date",{
  expect_true(is.na(pin_to_date(as.pin(c(NA,"198501169885")))[1]))
  expect_false(is.na(pin_to_date(as.pin(c(NA,"198501169885")))[2]))
  suppressWarnings(expect_true(all(is.na(pin_to_date(pin = c("19640823C234", "18640883D224"))))))
})
