context("Single data frame")

test_that("can read folder and store as multiple data frame", {
  folder_1 <- importAllSdf("MultipleFileTypes")
  folder_2 <- importAllSdf("MultipleFileTypes2")

  expect_equal(folder_1, folder_2)
})
