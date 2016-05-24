context("Test the threewords package")

# Set up test environment. To set the key:
# Sys.setenv("threewords_key", key_goes_here)
key <- Sys.getenv("threewords_key")

test_that("Individual word trios can be resolved to positions", {
  skip_on_cran()
  result <- from_words(key = key, words = c("turnip","basil","fruit"))
  expect_true(is.list(result))
  expect_equal(names(result), c("type","words","position","language"))
  expect_true(is.numeric(result$position))
})

test_that("Individual positions can be resolved to word trios", {
  skip_on_cran()
  result <- from_position(key = key, positions = c(6.385336, -36.293769))
  expect_true(is.list(result))
  expect_equal(names(result), c("words","position","language"))
  expect_true(is.numeric(result$position))
})


test_that("Multiple word trios can be resolved to positions", {
  skip_on_cran()
  result <- from_words(key = key, words = list(c("turnip","basil","fruit"), c("turnip","basil","fruit")))
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_equal(names(result[[1]]), c("type","words","position","language"))
  expect_true(is.numeric(result[[1]]$position))
})

test_that("Multiple coordinate sets can be resolved to positions", {
  skip_on_cran()
  result <- from_position(key = key, positions = list(c(6.385336, -36.293769), c(6.385336, -36.293769)))
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_equal(names(result[[1]]), c("words","position","language"))
  expect_true(is.numeric(result[[1]]$position))
})

test_that("Errors are correctly thrown",{
  expect_error(from_position(key = "ASDMASKMASLKMASLDMASLKMDASL", positions = c(6.385336, -36.293769)),
               regeexp = "missing or invalid key")
})