context("Test minicanada()")

test_that("Canada is drawn correctly", {
  expect_silent(minicanada(canada_abb, 1:13))
  expect_silent(minicanada(canada_abb, 1:13, border_colors = rep("red", 13), 
                        pt_names = FALSE))
  expect_silent(minicanada(canada_abb, rep("black", 13), border_colors = rep("white", 13), 
                        pt_name_colors = rep("yellow", 13), pt_name_cex = .5, 
                        font = "serif"))
})

test_that("minicanada() throws appropriate errors", {
  expect_error(minicanada(c(canada_abb[-1], "NU"), 1:13))
  expect_error(minicanada(canada_abb[1:10], 1:13))
  expect_error(minicanada(canada_abb, 1:10))
  expect_error(minicanada(canada_abb, 1:13, border_colors = 1:10))
  expect_error(minicanada(canada_abb, 1:13, pt_names = NA))
  expect_error(minicanada(canada_abb, 1:13, pt_name_colors = 1:10))
})