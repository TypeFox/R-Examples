context("oin S3 methods")

x <- as.oin(c("556000-4615", "232100-0156", "802002-4280"))

test_that("class attribut preserved", {
  expect_is(x[1], "oin")
  expect_is(`[.oin`(x, 1), "oin")
  expect_is(rep.oin(x, 3), "oin")
  expect_is({x[1] <- "232100-0156"; x}, "oin")
  expect_true(nchar(x[1])==11)
  expect_warning(x[1] <- "hejsan svejsan")
  expect_true(is.na(x[1]))
  expect_output(print(x), regexp = "Organizational identity number")
  expect_output(print(x), regexp = "556000-4615|232100-0156|802002-4280")
  expect_is(as.data.frame(x), "data.frame")
  expect_is(as.data.frame(x)[[1]], "oin")
  })
