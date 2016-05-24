context("Tilde")

test_that("can read folder and store as multiple data frame", {
  folder_1 <- importTilde("Tilde")
  folder_2 <- importTilde("Tilde2")

  expect_equal(folder_1, folder_2)
})
