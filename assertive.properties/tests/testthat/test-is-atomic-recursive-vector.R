test_that("test.is_atomic.array.returns_true", {
  expect_true(is_atomic(array()))
})

test_that("test.is_atomic.complex.returns_true", {
  expect_true(is_atomic(complex()))
})

test_that("test.is_atomic.integer.returns_true", {
  expect_true(is_atomic(integer()))
})

test_that("test.is_atomic.logical.returns_true", {
  expect_true(is_atomic(logical()))
})

test_that("test.is_atomic.matrix.returns_true", {
  expect_true(is_atomic(matrix()))
})

test_that("test.is_atomic.null.returns_true", {
  expect_true(is_atomic(NULL))
})

test_that("test.is_atomic.numeric.returns_true", {
  expect_true(is_atomic(numeric()))
})

test_that("test.is_atomic.raw.returns_true", {
  expect_true(is_atomic(raw()))
})

test_that("test.is_atomic.something_recursive.returns_false", {
  expect_false(is_atomic(list()))
})


test_that("test.is_nested.atomic.returns_false", {
  expect_false(is_nested(1:5))
})

test_that("test.is_nested.list_of_vectors.returns_false", {
  expect_false(is_nested(list(1:5)))
})

test_that("test.is_nested.list_of_lists.returns_true", {
  expect_true(is_nested(list(list(1:5))))
})


test_that("test.is_non_nested.atomic.returns_true", {
  expect_true(is_non_nested(1:5))
})

test_that("test.is_non_nested.list_of_vectors.returns_true", {
  expect_true(is_non_nested(list(1:5)))
})

test_that("test.is_non_nested.list_of_lists.returns_false", {
  expect_false(is_non_nested(list(list(1:5))))
})


test_that("test.is_recursive.a_call.returns_true", {
  expect_true(is_recursive(call("sin", "pi")))
})

test_that("test.is_recursive.a_data.frame.returns_true", {
  expect_true(is_recursive(data.frame()))
})

test_that("test.is_recursive.a_formula.returns_true", {
  expect_true(is_recursive(y ~ x))
})

test_that("test.is_recursive.a_function.returns_true", {
  expect_true(is_recursive(function() {
  }))
})

test_that("test.is_recursive.a_list.returns_true", {
  expect_true(is_recursive(list()))
})

test_that("test.is_recursive.an_expression.returns_true", {
  expect_true(is_recursive(expression()))
})

test_that("test.is_recursive.something_atomic.returns_false", {
  expect_false(is_recursive(1:10))
})

test_that("test.is_vector.character.returns_true", {
  expect_true(is_vector(character()))
})

test_that("test.is_vector.complex.returns_true", {
  expect_true(is_vector(complex()))
})

test_that("test.is_vector.expression.returns_true", {
  expect_true(is_vector(expression()))
})

test_that("test.is_vector.integer.returns_true", {
  expect_true(is_vector(integer()))
})

test_that("test.is_vector.list.returns_true", {
  expect_true(is_vector(list()))
})

test_that("test.is_vector.logical.returns_true", {
  expect_true(is_vector(logical()))
})

test_that("test.is_vector.not_a_vector.returns_false", {
  expect_false(is_vector(matrix()))
})

test_that("test.is_vector.numeric.returns_true", {
  expect_true(is_vector(numeric()))
})

test_that("test.is_vector.raw.returns_true", {
  expect_true(is_vector(raw()))
}) 
