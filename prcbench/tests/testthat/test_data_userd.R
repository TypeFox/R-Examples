context("Data: User-defeinded dataset")
# Test create_usrdata
#

test_that("create_usrdata: bench - test_type", {
  expect_error(create_usrdata("b", c(0.1, 0.2), c(1, 0)), NA)
  expect_error(create_usrdata("ben", c(0.1, 0.2), c(1, 0)), NA)

  expect_error(create_usrdata("bena", c(0.1, 0.2), c(1, 0)),
               "Invalid test_type")
})

test_that("create_usrdata: bench - scores", {
  expect_error(create_usrdata("bench", scores = c("0.1", "0.2"),
                              labels = c(1, 0)), "scores is not a numeric")

  expect_error(create_usrdata("bench", scores = 0.1,
                              labels = 1), "not greater than 1")
})

test_that("create_usrdata: bench - labels", {
  expect_error(create_usrdata("bench", scores = c(0.1, 0.2),
                              labels = c("1", "0")), "labels is not a numeric")

  expect_error(create_usrdata("bench", scores = c(0.1, 0.2),
                              labels = 1), "not greater than 1")

  expect_error(create_usrdata("bench", scores = c(0.1, 0.2),
                              labels = c(1, 1)), "not equal to 2")
})

test_that("create_usrdata: bench - scores and labels", {
  testset <- create_usrdata("bench", scores = c(0.1, 0.2),
                            labels = c(1, 0))[[1]]

  expect_true(is(testset, "TestDataB"))
  expect_true(is(testset, "R6"))

  expect_equal(testset$get_scores(), c(0.1, 0.2))
  expect_equal(testset$get_labels(), c(1, 0))

  expect_error(create_usrdata("bench", scores = c(0.1, 0.2),
                              labels = c(1, 0, 0)), "not equal to length")
})

test_that("create_usrdata: bench - tsname", {
  testset1 <- create_usrdata("bench", scores = c(0.1, 0.2),
                             labels = c(1, 0))[[1]]
  expect_equal(testset1$get_tsname(), "usr")

  testset2 <- create_usrdata("bench", scores = c(0.1, 0.2), labels = c(1, 0),
                             tsname = "m1")[[1]]
  expect_equal(testset2$get_tsname(), "m1")
})

test_that("create_usrdata: curve - test_type", {
  expect_error(create_usrdata("c", c(0.1, 0.2), c(1, 0), "A",
                              c(0, 1.0), c(0, 0.5)), NA)
  expect_error(create_usrdata("cur", c(0.1, 0.2), c(1, 0), "A",
                              c(0, 1.0), c(0, 0.5)), NA)

  expect_error(create_usrdata("cure", c(0.1, 0.2), c(1, 0), "A",
                              c(0, 1.0), c(0, 0.5)), "Invalid test_type")
})

test_that("create_usrdata: curve - scores", {
  expect_error(create_usrdata("curve", scores = c("0.1", "0.2"),
                              labels = c(1, 0)), "scores is not a numeric")

  expect_error(create_usrdata("curve", scores = 0.1,
                              labels = 1), "not greater than 1")
})

test_that("create_usrdata: curve - labels", {
  expect_error(create_usrdata("curve", scores = c(0.1, 0.2),
                              labels = c("1", "0")), "labels is not a numeric")

  expect_error(create_usrdata("curve", scores = c(0.1, 0.2),
                              labels = 1), "not greater than 1")

  expect_error(create_usrdata("curve", scores = c(0.1, 0.2),
                              labels = c(1, 1)), "not equal to 2")
})

test_that("create_usrdata: curve - scores and labels", {
  testset <- create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                            "A", c(0, 1.0), c(0, 0.5))[[1]]

  expect_true(is(testset, "TestDataB"))
  expect_true(is(testset, "R6"))

  expect_equal(testset$get_scores(), c(0.1, 0.2))
  expect_equal(testset$get_labels(), c(1, 0))

  expect_error(create_usrdata("curve", scores = c(0.1, 0.2),
                              labels = c(1, 0, 0), "A", c(0, 1.0), c(0, 0.5)),
               "not equal to length")
})

test_that("create_usrdata: curve - tsname", {
  testset1 <- create_usrdata("curve", scores = c(0.1, 0.2),
                             labels = c(1, 0), NULL, c(0, 1.0), c(0, 0.5))[[1]]
  expect_equal(testset1$get_tsname(), "usr")

  testset2 <- create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                             tsname = "m1", c(0, 1.0), c(0, 0.5))[[1]]
  expect_equal(testset2$get_tsname(), "m1")
})

test_that("create_usrdata: curve - base points", {
  testset <- create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                            base_x = c(0.13, 0.2), base_y = c(0.5, 0.6))[[1]]

  expect_true(is(testset, "TestDataC"))
  expect_true(is(testset, "TestDataB"))
  expect_true(is(testset, "R6"))

  expect_equal(testset$get_scores(), c(0.1, 0.2))
  expect_equal(testset$get_labels(), c(1, 0))
  expect_equal(testset$get_basepoints_x(), c(0.13, 0.2))
  expect_equal(testset$get_basepoints_y(), c(0.5, 0.6))

  expect_error(create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                              base_x = c(0.13, 0.2), base_y = c(0.5, 0.6, 1)),
               "not equal to length")
})

test_that("create_usrdata: curve - base_x", {
  expect_error(create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                              base_x = c("0.13", "0.2"), base_y = c(0.5, 0.6)),
               "base_x is not a numeric")

  expect_error(create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                              base_x = c(0.13, 1.1), base_y = c(0.5, 0.6)),
               "base_x <= 1 are not true")

  expect_error(create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                              base_x = c(-0.13, 1.1), base_y = c(0.5, 0.6)),
               "base_x >= 0 are not true")
})

test_that("create_usrdata: curve - base_y", {
  expect_error(create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                              base_x = c(0.13, 0.2), base_y = c("0.5", "0.6")),
               "base_y is not a numeric")

  expect_error(create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                              base_x = c(0.13, 0.2), base_y = c(1.5, 0.6)),
               "base_y <= 1 are not true")

  expect_error(create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                              base_x = c(0.13, 0.2), base_y = c(-0.5, 0.6)),
               "base_y >= 0 are not true")
})

test_that("create_usrdata: curve - text position", {

  testset2 <- create_usrdata("curve", scores = c(0.1, 0.2), labels = c(1, 0),
                             base_x = c(0.13, 0.2), base_y = c(0.5, 0.6),
                             text_x = 0.75, text_y = 0.85)[[1]]

  expect_true(is(testset2, "TestDataC"))
  expect_true(is(testset2, "TestDataB"))
  expect_true(is(testset2, "R6"))

  expect_equal(testset2$get_scores(), c(0.1, 0.2))
  expect_equal(testset2$get_labels(), c(1, 0))
  expect_equal(testset2$get_basepoints_x(), c(0.13, 0.2))
  expect_equal(testset2$get_basepoints_y(), c(0.5, 0.6))
  expect_equal(testset2$get_textpos_x(), 0.75)
  expect_equal(testset2$get_textpos_y(), 0.85)
})
