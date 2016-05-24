test_that("test.is_relistable.a_relistable_object.returns_true", {
  expect_true(is_relistable(as.relistable(list(1, 2, 3))))
})

test_that("test.is_relistable.a_list.returns_false", {
  expect_false(is_relistable(list(1, 2, 3)))
})
