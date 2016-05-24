context("pin S3 methods")

x <- as.pin(c(191212121212, 198505043334))

test_that("class attribut preserved", {
  expect_is(x[1], "pin")
  expect_is(`[.pin`(x, 1), "pin")
  expect_is(rep.pin(x, 3), "pin")
  expect_equal(length(rep.pin(x, 3)), 6)
  expect_is({x[1] <- "8505043334"; x}, "pin")
  expect_true(nchar(x[1])==12)
  expect_warning(x[1] <- "hejsan svejsan")
  expect_true(is.na(x[1]))
  expect_output(print(x), regexp = "Personal identity number")
  expect_output(print(x), regexp = "191212121212|198505043334")
  expect_is(as.data.frame(x), "data.frame")
  expect_is(as.data.frame(x)[[1]], "pin")
  })
