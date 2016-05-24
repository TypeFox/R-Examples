test_that("test.is_an_empty_string.empty_character.returns_false", {
  expect_false(is_an_empty_string(character()))
})

test_that("test.is_an_empty_string.empty_string.returns_true", {
  expect_true(is_an_empty_string(""))
})

test_that("test.is_an_empty_string.non_empty_string.returns_false", {
  expect_false(is_an_empty_string("a"))
})

test_that("test.is_an_empty_string.na_character_.returns_false", {
  expect_false(is_an_empty_string(NA_character_))
})


test_that("test.is_a_non_empty_string.empty_character.returns_false", {
  expect_false(is_a_non_empty_string(character()))
})

test_that("test.is_a_non_empty_string.empty_string.returns_false", {
  expect_false(is_a_non_empty_string(""))
})

test_that("test.is_a_non_empty_string.non_empty_string.returns_true", {
  expect_true(is_a_non_empty_string("a"))
})

test_that("test.is_a_non_empty_string.na_character_.returns_true", {
  expect_true(is_a_non_empty_string(NA_character_))
})


test_that("test.is_a_missing_or_empty_string.empty_character.returns_false", {
  expect_false(is_a_missing_or_empty_string(character()))
})

test_that("test.is_a_missing_or_empty_string.empty_string.returns_true", {
  expect_true(is_a_missing_or_empty_string(""))
})

test_that("test.is_a_missing_or_empty_string.non_empty_string.returns_false", {
  expect_false(is_a_missing_or_empty_string("a"))
})

test_that("test.is_a_missing_or_empty_string.na_character_.returns_false", {
  expect_true(is_a_missing_or_empty_string(NA_character_))
})


test_that("test.is_a_non_missing_nor_empty_string.empty_character.returns_false", {
  expect_false(is_a_non_missing_nor_empty_string(character()))
})

test_that("test.is_a_non_missing_nor_empty_string.empty_string.returns_false", {
  expect_false(is_a_non_missing_nor_empty_string(""))
})

test_that("test.is_a_non_missing_nor_empty_string.non_empty_string.returns_true", {
  expect_true(is_a_non_missing_nor_empty_string("a"))
})

test_that("test.is_a_non_missing_nor_empty_string.na_character_.returns_true", {
  expect_false(is_a_non_missing_nor_empty_string(NA_character_))
})

