
context("suppressAll")

test_that("suppressAll", {
  expect_equal(suppressAll(123), 123)
  expect_true(prints_text("123")({print(123);0})$passed)
  expect_false(prints_text("123")(suppressAll({print(123);0}))$passed)
  
  #todo: do later when testhat is fixed
  expect_true(gives_warning("123")({warning(123);0})$passed)
  #expect_false(gives_warning("123")(suppressAll({warning(123);0}))$passed)

  #todo: do later when testhat is fixed
  expect_true(shows_message("123")({message(123);0})$passed)
  #expect_false(shows_message("123")(suppressAll({message(123);0}))$passed)
})

