context("Test in_condition")

test_that("in_condition works as expected", {
  expect_identical(
    in_condition("COL_1", 1:3), "COL_1  in (1, 2, 3)")
  expect_identical(
    in_condition("COL_1", 1:3, "not"), "COL_1 not in (1, 2, 3)")
  expect_identical(
    in_condition("COL_1", LETTERS[2:3]), "COL_1  in ('B', 'C')")
  expect_identical(
    in_condition("COL_1", LETTERS[2:3], "not"), "COL_1 not in ('B', 'C')")
  expect_identical(
    in_condition("COL_1", ""), "COL_1  in ('')")
})

test_that("in_condition checks choices", {
  expect_error(in_condition("COL_1", as.POSIXct(Sys.Date())))
  expect_error(in_condition("COL_1", as.Date(NA)))
  expect_error(in_condition("COL_1", NA_character_))
  expect_error(in_condition("COL_1", integer(0)))
  expect_error(in_condition("COL_1", NULL))
  expect_error(in_condition("COL_1", 0.1))
  expect_error(in_condition("COL_1", "don't do this"))
  expect_error(in_condition("COL_1", 'or "this"'))
  expect_error(in_condition("COL_1", "or \"this"))
  expect_error(in_condition("COL_1", "or \'this"))
})

test_that("in_condition checks column_names", {
  expect_error(in_condition("", 1:3))
  expect_error(in_condition(NA, 1:3))
  expect_error(in_condition(1L, 1:3))
  expect_error(in_condition("'wrong'", 1:3))
  expect_error(in_condition("wrong wrong", 1:3))
  expect_error(in_condition("123", 1:3))
  expect_error(in_condition("ABC$", 1:3))
})


