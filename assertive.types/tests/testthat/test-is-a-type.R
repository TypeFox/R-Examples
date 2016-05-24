test_that("test.is_a_bool.a_vector.returns_false", {
  x <- c(TRUE, FALSE)
  actual <- is_a_bool(x)
  expect_false(actual)
  expect_equal(
    cause(actual), 
    noquote("x has length 2, not 1.")
  )
})

test_that("test.is_a_bool.empty_logical.returns_false", {
  x <- logical()
  actual <- is_a_bool(x)
  expect_false(actual)
  expect_equal(
    cause(actual), 
    noquote("x has length 0, not 1.")
  )
})

test_that("test.is_a_bool.false.returns_true", {
  expect_true(is_a_bool(FALSE))
})

test_that("test.is_a_bool.na.returns_true", {
  expect_true(is_a_bool(NA))
})

test_that("test.is_a_bool.true.returns_true", {
  expect_true(is_a_bool(TRUE))
})

test_that("test.is_a_complex.1.returns_false", {
  x <- 1L
  actual <- is_a_complex(x)
  expect_false(actual)
  expect_equal(
    cause(actual), 
    noquote("x is not of type 'complex'; it has class 'integer'.")
  )
})

test_that("test.is_a_complex.1_plus_0i.returns_true", {
  expect_true(is_a_complex(1 + (0 + (0 + (0 + (0+0i))))))
})

test_that("test.is_a_complex.1i.returns_true", {
  expect_true(is_a_complex(0 + (0 + (0 + (0+1i)))))
})

test_that("test.is_a_complex.a_vector.returns_false", {
  x <- c(0 + (0 + (0 + (0+1i))), 0 + (0 + (0 + (0+2i))))
  actual <- is_a_complex(x)
  expect_false(actual)
  expect_equal(
    cause(actual), 
    noquote("x has length 2, not 1.")
  )
})

test_that("test.is_a_complex.empty_complex.returns_false", {
  expect_false(is_a_complex(complex()))
})

test_that("test.is_a_complex.na_complex_.returns_true", {
  expect_true(is_a_complex(NA_complex_))
})

test_that("test.is_a_number.1.returns_true", {
  expect_true(is_a_number(1))
})

test_that("test.is_a_number.1L.returns_true", {
  expect_true(is_a_number(1L))
})

test_that("test.is_a_number.a_vector.returns_false", {
  expect_false(is_a_number(1:10))
})

test_that("test.is_a_number.empty_numeric.returns_false", {
  expect_false(is_a_number(numeric()))
})

test_that("test.is_a_number.Inf.returns_true", {
  expect_true(is_a_number(Inf))
})

test_that("test.is_a_number.na_real_.returns_true", {
  expect_true(is_a_number(NA_real_))
})

test_that("test.is_a_raw.a_raw.returns_true", {
  expect_true(is_a_raw(as.raw(1)))
})

test_that("test.is_a_raw.a_vector.returns_false", {
  expect_false(is_a_raw(as.raw(1:10)))
})

test_that("test.is_a_raw.empty_raw.returns_false", {
  expect_false(is_a_raw(raw()))
})

test_that("test.is_a_string.a_vector.returns_false", {
  expect_false(is_a_string(c("foo", "bar")))
})

test_that("test.is_a_string.empty_character.returns_false", {
  expect_false(is_a_string(character()))
})

test_that("test.is_a_string.empty_string.returns_true", {
  expect_true(is_a_string(""))
})

test_that("test.is_a_string.foo.returns_true", {
  expect_true(is_a_string("foo"))
})

test_that("test.is_a_string.na.returns_true", {
  expect_true(is_a_string(NA_character_))
})

test_that("test.is_an_integer.1L.returns_true", {
  expect_true(is_an_integer(1L))
})

test_that("test.is_an_integer.a_vector.returns_false", {
  expect_false(is_an_integer(1L:2L))
})

test_that("test.is_an_integer.empty_integer.returns_false", {
  expect_false(is_an_integer(integer()))
})

test_that("test.is_an_integer.floating_point.returns_false", {
  expect_false(is_an_integer(1))
})

test_that("test.is_an_integer.minus_1L.returns_true", {
  expect_true(is_an_integer(-1L))
})

test_that("test.is_an_integer.na.returns_true", {
  expect_true(is_an_integer(NA_integer_))
}) 

test_that("test.is_inherited_from.x_inherited_from_lowest_class.returns_true",  {
  x <- structure(1:5, class = c("foo", "bar"))
  expect_true(is_inherited_from(x, "foo"))   
})

test_that("test.is_inherited_from.x_inherited_from_highest_class.returns_true",  {
  x <- structure(1:5, class = c("foo", "bar"))
  expect_true(is_inherited_from(x, "bar"))   
})

test_that("test.is_inherited_from.x_only_inherited_from_some_classes.returns_true",  {
  x <- structure(1:5, class = c("foo", "bar"))
  expect_true(is_inherited_from(x, c("foo", "baz")))   
})

test_that("test.is_inherited_from.x_not_inherited",  {
  x <- structure(1:5, class = c("foo", "bar"))
  expect_false(actual <- is_inherited_from(x, "baz"))
  expect_equal(
    cause(actual),
    noquote("x does not inherit from the class baz. It has class foo, bar.")
  )
})


test_that("test.is_inherited_from.x_not_inherited_from_multiple",  {
  x <- structure(1:5, class = c("foo", "bar"))
  expect_false(actual <- is_inherited_from(x, c("baz", "quux")))
  expect_equal(
    cause(actual),
    noquote("x does not inherit from any of the classes baz, quux. It has class foo, bar.")
  )
})

