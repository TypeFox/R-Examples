context("XLSX")

test_that("can read folder and store as multiple data frame", {
  folder_1 <- importExcel("XLSX")
  folder_2 <- importExcel("XLSX2")

  expect_equal(folder_1, folder_2)
})
