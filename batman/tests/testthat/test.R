context("Test logical conversion functionality")

test_that("Simple logicals (including NAs) can be handled", {
  expect_that(to_logical(c("true","t","y","none","flkargfs")), equals(c(TRUE, TRUE, TRUE, FALSE, NA)))
})

test_that("Custom logical can be handled", {
  expect_that(to_logical(c("true","t","y","none","flkargfs","blargh"),
                         custom_true = "blargh"), equals(c(TRUE, TRUE, TRUE, FALSE, NA, TRUE)))
})

test_that("Custom logical can be handled", {
  expect_that(to_logical(c("true","t","y","none","flkargfs","blargh"),
                        custom_false = "blargh"), equals(c(TRUE, TRUE, TRUE, FALSE, NA, FALSE)))
})

test_that("Language codes can be retrieved", {
  
  result <- get_languages()
  expect_that(is.vector(result, mode = "character"), equals(TRUE))
  expect_that("en" %in% result, equals(TRUE))
})

test_that("Logicals in non-english languages can be handled", {
  values <- c("true","blargh","flargh", NA, "不對")
  expect_that(to_logical(values, language = "zh"), equals(c(NA,NA,NA,NA,FALSE)))
})
