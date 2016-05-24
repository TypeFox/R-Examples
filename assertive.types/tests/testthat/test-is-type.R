test_that("test.is_array.a_data.frame.returns_false", {
  expect_false(is_array(data.frame(x = 1:5)))
})

test_that("test.is_array.a_matrix.returns_true", {
  expect_true(is_array(matrix()))
})

test_that("test.is_array.a_vector.returns_false", {
  expect_false(is_array(1:10))
})

test_that("test.is_array.an_array.returns_true", {
  expect_true(is_array(array()))
})

test_that("test.is_atomic.not_a_call.returns_false", {
  expect_false(is_call(expression(sin(pi))))
})

test_that("test.is_call.a_call.returns_true", {
  expect_true(is_call(call("sin", "pi")))
})

test_that("test.is_character.character_vector.returns_true", {
  expect_true(is_character(letters))
})

test_that("test.is_character.NA_character_.returns_true", {
  expect_true(is_character(NA_character_))
})

test_that("test.is_character.not_a_character_vector.returns_false", {
  expect_false(is_character(1:10))
})

test_that("test.is_complex.1.returns_false", {
  expect_false(is_complex(1L))
})

test_that("test.is_complex.1_plus_0i.returns_true", {
  expect_true(is_complex(1 + (0 + (0 + (0 + (0+0i))))))
})

test_that("test.is_complex.1i.returns_true", {
  expect_true(is_complex(0 + (0 + (0 + (0+1i)))))
})

test_that("test.is_complex.na_complex_.returns_true", {
  expect_true(is_complex(NA_complex_))
})

test_that("test.is_data.frame.a_data.frame.returns_true", {
  expect_true(is_data.frame(data.frame(x = 1:5)))
})

test_that("test.is_data.frame.not_a_data.frame.returns_false", {
  expect_false(is_data.frame(list(x = 1:5)))
})

test_that("test.is_environment.an_environment.returns_true", {
  expect_true(is_environment(new.env()))
})

test_that("test.is_environment.base_environment.returns_true", {
  expect_true(is_environment(baseenv()))
})

test_that("test.is_environment.global_environment.returns_true", {
  expect_true(is_environment(globalenv()))
})

test_that("test.is_environment.not_an_environment.returns_false", {
  expect_false(is_environment(list()))
})

test_that("test.is_expression.an_environment.returns_true", {
  expect_true(is_expression(expression(sin(pi))))
})

test_that("test.is_expression.not_an_expression.returns_false", {
  expect_false(is_expression(call("sin", "pi")))
})

test_that("test.is_factor.a_factor.returns_true", {
  expect_true(is_factor(factor(letters)))
})

test_that("test.is_factor.an_ordered_factor.returns_true", {
  expect_true(is_factor(ordered(letters, levels = letters)))
})

test_that("test.is_factor.not_a_factor.returns_false", {
  expect_false(is_factor(letters))
})

test_that("test.is_function.a_function.returns_true", {
  expect_true(is_function(function() {
  }))
})

test_that("test.is_function.a_primitive_function.returns_true", {
  expect_true(is_function(sqrt))
})

test_that("test.is_function.not_a_function.returns_false", {
  expect_false(is_function(call("sin", "pi")))
})

test_that("test.is_integer.an_integer_vector.returns_true", {
  expect_true(is_integer(1L:10L))
})

test_that("test.is_integer.na_integer_.returns_true", {
  expect_true(is_integer(NA_integer_))
})

test_that("test.is_integer.not_an_integer.returns_false", {
  expect_false(is_integer(pi:10))
})

test_that("test.is_language.a_call.returns_true", {
  expect_true(is_language(call("sin", "pi")))
})

test_that("test.is_language.a_name.returns_true", {
  expect_true(is_language(as.name("foo")))
})

test_that("test.is_language.an_expression.returns_true", {
  expect_true(is_language(expression(sin(pi))))
})

test_that("test.is_language.not_a_language_object.returns_false", {
  expect_false(is_language(sin))
})

test_that("test.is_list.a_list.returns_true", {
  expect_true(is_list(list(1, 2, 3)))
})

test_that("test.is_list.an_atomic_vector.returns_false", {
  expect_false(is_list(1:10))
})

test_that("test.is_list.null.returns_false", {
  expect_false(is_list(NULL))
})

test_that("test.is_logical.a_logical_vector.returns_true", {
  expect_true(is_logical(c(TRUE, FALSE)))
})

test_that("test.is_logical.na.returns_true", {
  expect_true(is_logical(NA))
})

test_that("test.is_logical.not_a_logical.returns_false", {
  expect_false(is_logical(1:10))
})

test_that("test.is_matrix.a_data.frame.returns_false", {
  expect_false(is_matrix(data.frame(x = 1:5)))
})

test_that("test.is_matrix.a_matrix.returns_true", {
  expect_true(is_matrix(matrix()))
})

test_that("test.is_matrix.a_vector.returns_false", {
  expect_false(is_matrix(1:10))
})

test_that("test.is_matrix.an_array.returns_false", {
  expect_false(is_matrix(array()))
})

test_that("test.is_name.a_name.returns_true", {
  expect_true(is_name(as.name("foo")))
})

test_that("test.is_name.not_a_name.returns_false", {
  expect_false(is_name(call("sin", "pi")))
})

test_that("test.is_numeric.a_numeric_vector.returns_true", {
  expect_true(is_numeric(1:10))
})

test_that("test.is_numeric.an_integer_vector.returns_true", {
  expect_true(is_numeric(1L:10L))
})

test_that("test.is_numeric.not_numeric.returns_false", {
  expect_false(is_numeric(c(TRUE, FALSE)))
})

test_that("test.is_ordered.an_ordered_factor.returns_true", {
  expect_true(is_ordered(ordered(letters, levels = letters)))
})

test_that("test.is_ordered.an_unordered_factor.returns_false", {
  expect_false(is_ordered(factor(letters)))
})

test_that("test.is_ordered.not_a_factor.returns_false", {
  expect_false(is_ordered(letters))
})

test_that("test.is_primitive.a_primitive_function.returns_true", {
  expect_true(is_primitive(sqrt))
})

test_that("test.is_primitive.a_regular_function.returns_false", {
  expect_false(is_primitive(function() {
  }))
})

test_that("test.is_primitive.not_a_function.returns_false", {
  expect_false(is_primitive(call("sin", "pi")))
})

test_that("test.is_raw.a_raw_vector.returns_true", {
  expect_true(is_raw(as.raw(1:10)))
})

test_that("test.is_raw.not_raw.returns_false", {
  expect_false(is_raw(c(TRUE, FALSE)))
})

test_that("test.is_s4.an_S4_instance.returns_true", {
  x <- getClass("MethodDefinition")
  expect_true(is_s4(x))
})

test_that("test.is_s4.not_an_S4_instance.returns_true", {
  expect_false(is_s4(1:10))
})

test_that("test.is_table.a_table.returns_true", {
  x <- table(sample(letters, 100, replace = TRUE))
  expect_true(is_table(x))
})

test_that("test.is_table.not_a_table.returns_false", {
  expect_false(is_table(1:10))
})

test_that("test.is2.1_to_5_is_list.returns_false", {
  expect_false(is2(1:5, "list"))
})

test_that("test.is2.1_to_5_is_nonsense.returns_false", {
  expect_false(is2(1:5, "a b c"))
})

test_that("test.is2.1_to_5_is_numeric.returns_true", {
  expect_true(is2(1:5, "numeric"))
}) 
