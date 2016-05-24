context("Tab")

test_that("can read folder and store as multiple data frame", {
  folder_1 <- importTab("Tab")
  folder_2 <- importTab("Tab2")

  expect_equal(folder_1, folder_2)
})
