context("Data: TestData")
# Test TestDataB
#      TestDataC
#

test_that("TestDataB - R6ClassGenerator", {
  expect_true(is(TestDataB, "R6ClassGenerator"))
  expect_equal(attr(TestDataB, "name"), "TestDataB_generator")

  expect_true(is.function(TestDataB$public_methods$get_tsname))
  expect_true(is.function(TestDataB$public_methods$get_scores))
  expect_true(is.function(TestDataB$public_methods$get_labels))
  expect_true(is.function(TestDataB$public_methods$get_fg))
  expect_true(is.function(TestDataB$public_methods$get_bg))
  expect_true(is.function(TestDataB$public_methods$get_fname))
  expect_true(is.function(TestDataB$public_methods$del_file))
})

test_that("TestDataB - R6", {
  data_obj <- TestDataB$new(c(0.1, 0.2, 0.3), c(0, 1, 1))

  expect_true(is(data_obj, "TestDataB"))
  expect_true(is(data_obj, "R6"))

  expect_true(is.function(data_obj[["get_tsname"]]))
  expect_true(is.function(data_obj[["get_scores"]]))
  expect_true(is.function(data_obj[["get_labels"]]))
  expect_true(is.function(data_obj[["get_fg"]]))
  expect_true(is.function(data_obj[["get_bg"]]))
  expect_true(is.function(data_obj[["get_fname"]]))
  expect_true(is.function(data_obj[["del_file"]]))
})

test_that("TestDataB - get_datname", {
  data_obj1 <- TestDataB$new(c(0.1, 0.2, 0.3), c(0, 1, 1))
  expect_true(is.na(data_obj1$get_tsname()))

  data_obj2 <- TestDataB$new(c(0.1, 0.2, 0.3), c(0, 1, 1), "m1")
  expect_equal(data_obj2$get_tsname(), "m1")
})

test_that("TestDataB - get_scores", {
  data_obj <- TestDataB$new(c(0.1, 0.2, 0.3), c(0, 1, 1))
  expect_equal(data_obj$get_scores(), c(0.1, 0.2, 0.3))
})

test_that("TestDataB - get_labels", {
  data_obj <- TestDataB$new(c(0.1, 0.2, 0.3), c(0, 1, 1))
  expect_equal(data_obj$get_labels(), c(0, 1, 1))
})

test_that("TestDataB - get_fg", {
  data_obj <- TestDataB$new(c(0.1, 0.2, 0.3, 0.4), c(0, 1, 1, 0))
  expect_equal(data_obj$get_fg(), c(0.2, 0.3))
})

test_that("TestDataB - get_bg", {
  data_obj <- TestDataB$new(c(0.1, 0.2, 0.3, 0.4), c(0, 1, 1, 0))
  expect_equal(data_obj$get_bg(), c(0.1, 0.4))
})

test_that("TestDataB - get_fname", {
  data_obj <- TestDataB$new(c(0.1, 0.2, 0.3, 0.4), c(0, 1, 1, 0))
  expect_true(file.exists(data_obj$get_fname()))
})

test_that("TestDataB - del_file", {
  data_obj <- TestDataB$new(c(0.1, 0.2, 0.3, 0.4), c(0, 1, 1, 0))
  fname <- data_obj$get_fname()

  data_obj$del_file()

  expect_true(!file.exists(fname))
  expect_true(is.na(data_obj$get_fname()))
})

test_that("TestDataC - R6ClassGenerator", {
  expect_true(is(TestDataC, "R6ClassGenerator"))
  expect_equal(attr(TestDataC, "name"), "TestDataC_generator")

  expect_true(is.function(TestDataC$public_methods$set_basepoints_x))
  expect_true(is.function(TestDataC$public_methods$set_basepoints_y))
  expect_true(is.function(TestDataC$public_methods$get_basepoints_x))
  expect_true(is.function(TestDataC$public_methods$get_basepoints_y))
  expect_true(is.function(TestDataC$public_methods$set_textpos_x))
  expect_true(is.function(TestDataC$public_methods$set_textpos_y))
  expect_true(is.function(TestDataC$public_methods$get_textpos_x))
  expect_true(is.function(TestDataC$public_methods$get_textpos_y))
})

test_that("TestDataC - R6", {
  data_obj <- TestDataC$new(c(0.1, 0.2, 0.3), c(0, 1, 1))

  expect_true(is(data_obj, "TestDataC"))
  expect_true(is(data_obj, "TestDataB"))
  expect_true(is(data_obj, "R6"))

  expect_true(is.function(data_obj[["set_basepoints_x"]]))
  expect_true(is.function(data_obj[["set_basepoints_y"]]))
  expect_true(is.function(data_obj[["get_basepoints_x"]]))
  expect_true(is.function(data_obj[["get_basepoints_y"]]))
  expect_true(is.function(data_obj[["set_textpos_x"]]))
  expect_true(is.function(data_obj[["set_textpos_y"]]))
  expect_true(is.function(data_obj[["get_textpos_x"]]))
  expect_true(is.function(data_obj[["get_textpos_y"]]))
})

test_that("TestDataC - get_datname", {
  data_obj1 <- TestDataC$new(c(0.1, 0.2, 0.3), c(0, 1, 1))
  expect_true(is.na(data_obj1$get_tsname()))

  data_obj2 <- TestDataC$new(c(0.1, 0.2, 0.3), c(0, 1, 1), "m1")
  expect_equal(data_obj2$get_tsname(), "m1")
})

test_that("TestDataC - get_scores", {
  data_obj <- TestDataC$new(c(0.1, 0.2, 0.3), c(0, 1, 1))
  expect_equal(data_obj$get_scores(), c(0.1, 0.2, 0.3))
})

test_that("TestDataC - basepoints", {
  data_obj <- TestDataC$new(c(0.1, 0.2, 0.3), c(0, 1, 1))
  data_obj$set_basepoints_x(c(0, 0.5, 1))
  data_obj$set_basepoints_y(c(0, 0.4, 1))

  expect_equal(data_obj$get_basepoints_x(), c(0, 0.5, 1))
  expect_equal(data_obj$get_basepoints_y(), c(0, 0.4, 1))
})

test_that("TestDataC - textpos", {
  data_obj <- TestDataC$new(c(0.1, 0.2, 0.3), c(0, 1, 1))
  data_obj$set_textpos_x(c(0.3, 0.4))
  data_obj$set_textpos_y(c(0.8, 0.9))

  expect_equal(data_obj$get_textpos_x(), c(0.3, 0.4))
  expect_equal(data_obj$get_textpos_y(), c(0.8, 0.9))
})
