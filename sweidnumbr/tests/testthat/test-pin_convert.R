
context("pin_convert")

test_that(desc="pin_convert",{
  expect_equal(pin_convert(pin = "198112189876", format = 1), expected = "198112189876")
  expect_equal(pin_convert(pin = "811218-9876", format = 3), expected = "198112189876")
  expect_equal(pin_convert(pin = "19811218-9876", format = 2), expected = "198112189876")
  expect_equal(pin_convert(pin = "8112189876", format = 4), expected = "198112189876")
  expect_equal(pin_convert(pin = "001218-0000", format = 3), expected = "200012180000")
  expect_is(pin_convert(pin = "001218-0000", format = 3), "character")
})
