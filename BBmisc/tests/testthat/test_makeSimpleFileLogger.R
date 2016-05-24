context("makeSimpleFileLogger")

test_that("makeSimpleFileLogger", {
  fn = tempfile()
  logger = makeSimpleFileLogger(fn)
  expect_identical(class(logger), "SimpleFileLogger")

  expect_equal(logger$getSize(), 0)
  msg1 = "xxx111xxx"
  logger$log(msg1)
  expect_true(grepl(msg1, readLines(fn)))
  expect_identical(msg1, logger$getMessages(1))
  expect_equal(logger$getSize(), 1)

  expect_true(file.exists(fn))
  logger$clear()
  expect_false(file.exists(fn))
})

test_that("message order",  {
  fn = tempfile()
  msg1 = "xxx111xxx"
  msg2 = "xxx222xxx"
  for (keep in c(0, 10)) {
    logger = makeSimpleFileLogger(fn)
    logger$log(msg1)
    logger$log(msg2)
    expect_identical(grepl("xxx[0-9]+xxx$", readLines(fn)), c(TRUE, TRUE))
    expect_identical(grepl("^xxx[0-9]+xxx$", logger$getMessages(2)), c(TRUE, TRUE))
    expect_identical(logger$getMessages(1), msg2)
    expect_identical(logger$getMessages(2), c(msg2, msg1))
  }
})
