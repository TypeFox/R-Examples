library(dplyr)

context("distance_join")

sepal_lengths <- data_frame(Sepal.Length = c(5, 6, 7), Sepal.Width = 1:3, Type = 1:3)

test_that("distance_inner_join works for Euclidean distance", {
  ret <- iris %>%
    distance_inner_join(sepal_lengths, max_dist = .25,
                        by = c("Sepal.Length", "Sepal.Width")) %>%
    mutate(distance = sqrt((Sepal.Length.x - Sepal.Length.y) ^ 2 +
                             (Sepal.Width.x - Sepal.Width.y) ^ 2))

  expect_true(nrow(ret) > 0)
  expect_true(all(ret$distance < .25))
})

test_that("distance_inner_join works for Manhattan distance", {
  ret2 <- iris %>%
    distance_inner_join(sepal_lengths, max_dist = .25,
                        by = c("Sepal.Length", "Sepal.Width"),
                        method = "manhattan") %>%
    mutate(distance = abs(Sepal.Length.x - Sepal.Length.y) +
                      abs(Sepal.Width.x - Sepal.Width.y))

  expect_true(nrow(ret2) > 0)
  expect_true(all(ret2$distance < .25))
})
