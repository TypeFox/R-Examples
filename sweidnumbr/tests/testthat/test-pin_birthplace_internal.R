
context("pin_birthplace_internal")


test_that(desc="birthplace",{
  expect_equal(pin_birthplace_internal(pin = "199000000100", birth_vector = as.character(1:100), birth_other_text = "Born after 31 december 1989"), expected = "Born after 31 december 1989")
  expect_equal(pin_birthplace_internal(pin = "194000000100", birth_vector = as.character(1:100), birth_other_text = "Born after 31 december 1989"), expected = "2")
  expect_equal(pin_birthplace_internal(pin = "194000009900", birth_vector = as.character(1:100), birth_other_text = "Born after 31 december 1989"), expected = "100")
  expect_is(pin_birthplace_internal(pin = "194000009900", birth_vector = as.character(1:100), birth_other_text = "Born after 31 december 1989"), "character")
})
