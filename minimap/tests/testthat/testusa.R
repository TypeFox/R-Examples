context("Test miniusa()")

test_that("The USA is drawn correctly", {
  expect_silent(miniusa(usa_abb, 1:51))
  expect_silent(miniusa(usa_abb, 1:51, border_colors = rep("red", 51), 
                        state_names = FALSE))
  expect_silent(miniusa(usa_abb, rep("black", 51), border_colors = rep("white", 51), 
                        state_name_colors = rep("yellow", 51), state_name_cex = .5, 
                        font = "serif"))
})

test_that("miniusa() throws appropriate errors", {
  expect_error(miniusa(c(usa_abb[-1], "DC"), 1:51))
  expect_error(miniusa(usa_abb[1:10], 1:51))
  expect_error(miniusa(usa_abb, 1:10))
  expect_error(miniusa(usa_abb, 1:51, border_colors = 1:10))
  expect_error(miniusa(usa_abb, 1:51, state_names = NA))
  expect_error(miniusa(usa_abb, 1:51, state_name_colors = 1:10))
})