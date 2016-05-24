context("Multi data frame")

test_that("can read folder and store as multiple data frame", {
  folder_1 <- importAllMdf("MultipleFileTypes")
  folder_2 <- importAllMdf("MultipleFileTypes2")

  expect_equal(folder_1, folder_2)
})
