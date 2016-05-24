
context("pin_coordn_correct")

test_that(desc="control number",{
  expect_equal(pin_coordn_correct(pin = "198112189876"), expected = "198112189876")
  expect_equal(pin_coordn_correct(pin = "198112789876"), expected = "198112189876")
  expect_equal(pin_coordn_correct(pin = "198112680000"), expected = "198112080000")
  expect_equal(pin_coordn_correct(pin = "198112880000"), expected = "198112280000")
  expect_is(pin_coordn_correct(pin = "198112880000"), "character")
})
