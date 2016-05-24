context("findSimilarName")

test_that("verifies findSimilarName", {
   x <- "xyz"
   dist1 <- "Ayz"
   dist2 <- "ABz"
   dist3 <- "ABC"
   actual <- DYM:::findSimilarName(x, c(dist1, dist2, dist3))
   # threshold defaults to 2
   expect_equal(actual, c(dist1, dist2))
})

# threshold ---------------------------------------------------------------

test_that("verifies findSimilarName with threshold = 3", {
   x <- "xyz"
   dist1 <- "Ayz"
   dist2 <- "ABz"
   dist3 <- "ABC"
   actual <- DYM:::findSimilarName(x, c(dist1, dist2, dist3), threshold=3)
   expect_equal(actual, c(dist1, dist2, dist3))
})

test_that("verifies findSimilarName with threshold = 1", {
   x <- "xyz"
   dist1 <- "Ayz"
   dist2 <- "ABz"
   dist3 <- "ABC"
   actual <- DYM:::findSimilarName(x, c(dist1, dist2, dist3), threshold=1)
   expect_equal(actual, dist1)
})

test_that("verifies findSimilarName with threshold = 0", {
   x <- "xyz"
   dist1 <- "Ayz"
   dist2 <- "ABz"
   dist3 <- "ABC"
   actual <- DYM:::findSimilarName(x, c(dist1, dist2, dist3), threshold=0)
   expect_equal(actual, character())
})

# NAs ---------------------------------------------------------------------

test_that("no result if x = NA", {
   actual <- DYM:::findSimilarName(NA, c("x", "y"))
   expect_equal(actual, character(0))
})

test_that("NA in names should be omitted", {
   actual <- DYM:::findSimilarName("N", c("x", NA, "y", "NA"))
   expect_equal(actual, c("x", "y", "NA"))
})

# max_n -------------------------------------------------------------------

test_that("verifies findSimilarName with max_n +ve and < n results", {
   x <- "tuvwxyz"
   # choices all distance 1 away from x, so all should match
   choices <- c("Auvwxyz", "tBvwxyz", "tuCwxyz", "tuvDxyz", "tuvwEyz", "tuvwxFz", "tuvwxyG")
   # -ve max_n means keep first max_n matches
   actual <- DYM:::findSimilarName(x, choices, max_n = 6)
   expect_equal(actual, choices[1:6])
})

test_that("verifies findSimilarName with max_n -ve and < n results", {
   x <- "tuvwxyz"
   # choices all distance 1 away from x, so all should match
   choices <- c("Auvwxyz", "tBvwxyz", "tuCwxyz", "tuvDxyz", "tuvwEyz", "tuvwxFz", "tuvwxyG")
   # -ve max_n means keep all but last max_n matches
   actual <- DYM:::findSimilarName(x, choices, max_n = -2)
   expect_equal(actual, choices[1:5])
})

# ignoreCase --------------------------------------------------------------

test_that("verifies findSimilarName with ignoreCase = FALSE", {
   x <- "xyz"
   changeCase <- c("Xyz", "XYz", "XYZ")
   changeLetter <- c("ayz", "abz", "abc")
   actual <- DYM:::findSimilarName(x, c(changeCase, changeLetter))
   expect_equal(actual, c(changeCase[1], changeLetter[1], changeCase[2], changeLetter[2]))
})

test_that("verifies findSimilarName with ignoreCase = TRUE", {
   x <- "xyz"
   changeCase <- c("Xyz", "XYz", "XYZ")
   changeLetter <- c("ayz", "abz", "abc")
   actual <- DYM:::findSimilarName(x, c(changeCase, changeLetter), ignoreCase = TRUE)
   expect_equal(actual, c(changeCase, changeLetter[1:2]))
})


