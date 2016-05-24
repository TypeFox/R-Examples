context("Print: S3 print")
# Test print.benchmark
#      print.evalcurve
#

test_that("print.benchmark", {
  testset <- create_testset("bench", "b10")
  toolset <- create_toolset(set_names = "crv5")
  res1<- run_benchmark(testset, toolset, times = 1)

  expect_error(capture.output(print(res1)), NA)

  expect_error(print(res1, digits = "2"), "digits is not a number")
})

test_that("print.evalcurve", {
  testset <- create_testset("curve", "c1")
  toolset <- create_toolset(set_names = "crv5")
  res1 <- run_evalcurve(testset, toolset)

  expect_error(capture.output(print(res1)), NA)
  expect_error(capture.output(print(res1, "summary")), NA)
  expect_error(capture.output(print(res1, "all")), NA)
  expect_error(capture.output(print(res1, "basepoints")), NA)
  expect_error(capture.output(print(res1, "predictions")), NA)
  expect_error(capture.output(print(res1, "rawsummary")), NA)

  expect_output(print(res1), "score")
  expect_output(print(res1, "summary"), "score")
  expect_output(print(res1, "all"), "testitem")
  expect_output(print(res1, "basepoints"), "testset")
  expect_output(print(res1, "predictions"), "x")
  expect_output(print(res1, "rawsummary"),  "lbl_pos_x")
})
