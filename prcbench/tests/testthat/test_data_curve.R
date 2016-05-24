context("Data: Testset for curve evaluation")
# Test create_testset
#

test_that("create_testset: test_type", {
  expect_error(create_testset("c", "c1"), NA)
  expect_error(create_testset("cur", "c1"), NA)

  expect_error(create_testset("cure", "c1"), "Invalid test_type")
})

test_that("create_testset: set_names", {
  expect_error(create_testset("curve", "c1"), NA)
  expect_error(create_testset("curve", "c2"), NA)
  expect_error(create_testset("curve", "c3"), NA)
  expect_error(create_testset("curve", "C1"), NA)
  expect_error(create_testset("curve", "C2"), NA)
  expect_error(create_testset("curve", "C3"), NA)
  expect_error(create_testset("curve", c("c1", "c2")), NA)

  expect_error(create_testset("curve", "c4"), "Invalid set_names")
  expect_error(create_testset("curve", "a1"), "Invalid set_names")
  expect_error(create_testset("curve", "abc"), "Invalid set_names")
  expect_error(create_testset("curve", c("c1", "a1")), "Invalid set_names")
})

test_that("create_testset: c1", {
  testset <- create_testset("curve", "c1")[[1]]

  expect_true(is(testset, "TestDataC"))
  expect_true(is(testset, "TestDataB"))
  expect_true(is(testset, "R6"))

  expect_equal(testset$get_tsname(), "c1")

  bp_x <- c(0, 0.25, 0.5, 0.75, 1, 1)
  bp_y <- c(1, 1, 1, 0.75, 0.6666666667, 0.5)
  expect_equal(testset$get_basepoints_x(), bp_x)
  expect_equal(testset$get_basepoints_y(), bp_y)

  expect_equal(testset$get_textpos_x(), 0.85)
  expect_equal(testset$get_textpos_y(), 0.9)
})

test_that("create_testset: c2", {
  testset <- create_testset("curve", "c2")[[1]]

  expect_true(is(testset, "TestDataC"))
  expect_true(is(testset, "TestDataB"))
  expect_true(is(testset, "R6"))

  expect_equal(testset$get_tsname(), "c2")

  bp_x <- c(0, 0.25, 0.5, 0.5, 0.75, 1)
  bp_y <- c(0.5, 0.5, 0.5, 0.3333333333, 0.4285714286, 0.5)
  expect_equal(testset$get_basepoints_x(), bp_x)
  expect_equal(testset$get_basepoints_y(), bp_y)

  expect_equal(testset$get_textpos_x(), 0.2)
  expect_equal(testset$get_textpos_y(), 0.65)
})

test_that("create_testset: c3", {
  testset <- create_testset("curve", "c3")[[1]]

  expect_true(is(testset, "TestDataC"))
  expect_true(is(testset, "TestDataB"))
  expect_true(is(testset, "R6"))

  expect_equal(testset$get_tsname(), "c3")

  bp_x <- c(0, 0, 0, 0.25, 0.5, 0.75, 1)
  bp_y <- c(0, 0, 0, 0.2, 0.3333333333, 0.4285714286, 0.5)
  expect_equal(testset$get_basepoints_x(), bp_x)
  expect_equal(testset$get_basepoints_y(), bp_y)

  expect_equal(testset$get_textpos_x(), 0.8)
  expect_equal(testset$get_textpos_y(), 0.2)
})

test_that("create_testset: c1, c2, c3", {
  testset <- create_testset("curve", c("c1", "c2", "c3"))

  expect_equal(testset[[1]]$get_tsname(), "c1")
  expect_equal(testset[[2]]$get_tsname(), "c2")
  expect_equal(testset[[3]]$get_tsname(), "c3")
})
