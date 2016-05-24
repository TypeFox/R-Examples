library(ggplot2)
library(dplyr)

context("stringdist_join")

# setup
d <- data_frame(cut2 = c("Idea", "Premiums", "Premiom",
                         "VeryGood", "VeryGood", "Faiir")) %>%
  mutate(type = row_number())

test_that("stringdist_inner_join works on a large df with multiples in each", {
  # create something with names close to the cut column in the diamonds dataset
  j <- stringdist_inner_join(diamonds, d, by = c(cut = "cut2"))

  result <- j %>%
    count(cut, cut2) %>%
    arrange(cut)

  expect_equal(as.character(result$cut), c("Fair", "Very Good", "Premium", "Premium", "Ideal"))
  expect_equal(result$cut2, c("Faiir", "VeryGood", "Premiom", "Premiums", "Idea"))

  expect_equal(sum(j$cut == "Premium"), sum(diamonds$cut == "Premium") * 2)
  expect_equal(sum(j$cut == "Very Good"), sum(diamonds$cut == "Very Good") * 2)
  expect_equal(sum(j$cut2 == "Premiom"), sum(diamonds$cut == "Premium"))

  vg <- j %>%
    filter(cut == "Very Good") %>%
    count(type)

  expect_equal(vg$type, c(4, 5))
  expect_equal(vg$n, rep(sum(diamonds$cut == "Very Good"), 2))

  expect_true(all(j$type[j$cut == "Faiir"] == 1))
})


d2 <- head(d, 3)
included <- c("Ideal", "Premium")
notin <- c("Fair", "Good", "Very Good")


test_that("stringdist_left_join works as expected", {
  result <- diamonds %>%
    stringdist_left_join(d2, by = c(cut = "cut2"))

  expect_true(all(is.na(result$cut2[result$cut %in% notin])))
  expect_equal(sum(result$cut %in% notin), sum(diamonds$cut %in% notin))

  expect_equal(sum(result$cut2 == "Premiom", na.rm = TRUE),
               sum(diamonds$cut == "Premium"))
  expect_equal(sum(result$cut2 == "Premiom", na.rm = TRUE),
               sum(result$cut2 == "Premiums", na.rm = TRUE))
})


d3 <- bind_rows(d2, data_frame(cut2 = "NewType", type = 4))

test_that("stringdist_right_join works as expected", {
  result <- diamonds %>%
    stringdist_right_join(d3, by = c(cut = "cut2"))

  expect_equal(sum(result$cut2 == "NewType"), 1)
  expect_equal(sum(is.na(result$cut)), 1)
  expect_true(all(is.na(result$cut[result$cut2 == "NewType"])))

  expect_equal(sum(result$cut2 == "Premiom", na.rm = TRUE),
               sum(diamonds$cut == "Premium"))
  expect_equal(sum(result$cut2 == "Premiom", na.rm = TRUE),
               sum(result$cut2 == "Premiums", na.rm = TRUE))
})


test_that("stringdist_full_join works as expected", {
  result <- diamonds %>%
    stringdist_full_join(d3, by = c(cut = "cut2"))

  expect_equal(sum(result$cut2 == "NewType", na.rm = TRUE), 1)
  expect_equal(sum(is.na(result$cut)), 1)
  expect_true(all(is.na(result$cut[result$cut2 == "NewType"])))

  expect_true(all(is.na(result$cut2[result$cut %in% notin])))
  expect_equal(sum(result$cut %in% notin), sum(diamonds$cut %in% notin))

  expect_equal(sum(result$cut2 == "Premiom", na.rm = TRUE),
               sum(diamonds$cut == "Premium"))
  expect_equal(sum(result$cut2 == "Premiom", na.rm = TRUE),
               sum(result$cut2 == "Premiums", na.rm = TRUE))
})



test_that("stringdist_semi_join works as expected", {
  result <- diamonds %>%
    stringdist_semi_join(d2, by = c(cut = "cut2"))

  expect_equal(sort(as.character(unique(result$cut))), included)

  expect_equal(nrow(result), sum(result$cut %in% included))

  expect_true(!("cut2" %in% colnames(result)))
})



test_that("stringdist_anti_join works as expected", {
  result <- diamonds %>%
    stringdist_anti_join(d2, by = c(cut = "cut2"))

  expect_equal(sort(as.character(unique(result$cut))), notin)

  expect_equal(nrow(result), sum(result$cut %in% notin))
})


test_that("stringdist_inner_join works with multiple match functions", {
  # setup
  d3 <- data_frame(cut2 = c("Idea", "Premiums", "Premiom",
                           "VeryGood", "VeryGood", "Faiir"),
                   carat2 = c(0, .5, 1, 1.5, 2, 2.5)) %>%
    mutate(type = row_number())

  sdist <- function(s1, s2) stringdist::stringdist(s1, s2) <= 1
  ndist <- function(n1, n2) abs(n1 - n2) < .25

  j <- diamonds %>%
    fuzzy_inner_join(d3, by = c(cut = "cut2", carat = "carat2"),
                     match_fun = list(sdist, ndist))

  result <- j %>%
    count(cut, cut2)

  expect_equal(as.character(result$cut), c("Fair", "Very Good", "Premium", "Premium", "Ideal"))
  expect_equal(result$cut2, c("Faiir", "VeryGood", "Premiom", "Premiums", "Idea"))

  expect_less_than(max(abs(j$carat - j$carat2)), .25)

  # give match_fun as a named list
  j_named <- diamonds %>%
    fuzzy_inner_join(d3, by = c(cut = "cut2", carat = "carat2"),
                     match_fun = list(carat = ndist, cut = sdist))

  expect_equal(j, j_named)
})


test_that("stringdist_join works with data frames without matches", {
  d <- data_frame(cut2 = c("Ideolll", "Premiumsss", "Premiomzzz",
                           "VeryVeryGood", "VeryVeryGood", "FaiirsFair")) %>%
    mutate(type = row_number())

  j1 <- stringdist_inner_join(diamonds, d, by = c(cut = "cut2"))
  expect_equal(nrow(j1), 0)
  expect_true(all(c("carat", "cut", "cut2", "type") %in% colnames(j1)))

  # check it works when column names are the same
  d2 <- rename(d, cut = cut2)
  j1_5 <- stringdist_inner_join(diamonds, d2, by = c(cut = "cut"))
  expect_equal(nrow(j1_5), 0)
  expect_true(all(c("carat", "cut.x", "cut.y", "type") %in% colnames(j1_5)))

  j2 <- stringdist_left_join(diamonds, d, by = c(cut = "cut2"))
  expect_equal(nrow(j2), nrow(diamonds))
  expect_true(all(is.na(j2$cut2)))

  j3 <- stringdist_right_join(diamonds, d, by = c(cut = "cut2"))
  expect_equal(nrow(j3), nrow(d))
  expect_true(all(is.na(j3$carat)))
  expect_true(all(is.na(j3$cut)))

  j4 <- stringdist_full_join(diamonds, d, by = c(cut = "cut2"))
  expect_equal(nrow(j4), nrow(diamonds), nrow(d))
  expect_true(all(is.na(j4$cut) | is.na(j4$cut2)))
  expect_true(all(is.na(j4$carat) | is.na(j4$type)))

  j5 <- stringdist_semi_join(diamonds, d, by = c(cut = "cut2"))
  expect_equal(nrow(j5), 0)
  expect_true(!("cut2" %in% colnames(diamonds)))
  expect_true(!("type" %in% colnames(diamonds)))
  expect_true("cut" %in% colnames(diamonds))
  expect_true("carat" %in% colnames(diamonds))

  j6 <- stringdist_anti_join(diamonds, d, by = c(cut = "cut2"))
  expect_equal(nrow(j6), nrow(diamonds))
  expect_true("cut" %in% colnames(diamonds))
  expect_true(!("cut2" %in% colnames(diamonds)))
})
