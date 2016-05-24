################################################################################
# Basic pattern analysis functions
################################################################################

context("basic pattern analysis functions")

test_that("get_pattern works properly", {
  x <- c(" 12-Ac$&abnd abn", "Male", "01/01/1999", 3.1415)
  pat1 <- get_pattern(x, show_ws = FALSE)
  pat2 <- get_pattern(x, show_ws = TRUE)
  expect_identical(pat1, c(" 99-Aa$&aaaa aaa", "Aaaa", "99/99/9999", "9.9999"))
  expect_identical(pat2, c("w99-Aa$&aaaawaaa", "Aaaa", "99/99/9999", "9.9999"))
})

test_that("basic_pattern_analysis and bpa work properly", {
  d <- data.frame(x = c(" 12-Ac$&abnd abn", "Male", "01/01/1999", 3.1415),
                  y = letters[1:4],
                  z = 101:104)

  # Data frame output
  pat <- basic_pattern_analysis(d)
  expect_is(pat, "data.frame")
  expect_equal(attributes(pat), attributes(d))

  # List output
  pat <- basic_pattern_analysis(d, unique_only = TRUE)
  expect_is(pat, "list")
  expect_equal(names(pat), names(d))

  expect_equal(bpa(d), basic_pattern_analysis(d))

})


test_that("match_pattern works properly", {
  phone <- c("123-456-7890", "456-7890", "123-4567", "456-7890")
  mp1 <- match_pattern(phone, pattern = "999-9999")
  mp2 <- match_pattern(phone, pattern = "999-9999", unique_only = TRUE)
  expect_equal(mp1, c("456-7890", "123-4567", "456-7890"))
  expect_equal(mp2, c("456-7890", "123-4567"))
})
