library(ggplot2)
library(dplyr)

context("difference_inner_join")

test_that("difference_inner_join works on a df with multiples in each", {
  sepal_lengths <- data_frame(Sepal.Length = c(5, 6, 7), Type = 1:3)

  j <- iris %>%
    difference_inner_join(sepal_lengths, max_dist = .5)

  expect_equal(min(j$Sepal.Length.x - j$Sepal.Length.y), -.5)
  expect_equal(max(j$Sepal.Length.x - j$Sepal.Length.y), .5)

  expect_true(all(j$Type[j$Sepal.Length.x < 5.5] == 1))
})
