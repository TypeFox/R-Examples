context("Test minimexico()")

test_that("Mexico is drawn correctly", {
  expect_silent(minimexico(mexico_abb, 1:32))
  expect_silent(minimexico(mexico_abb, 1:32, border_colors = rep("red", 32), 
                        estados_names = FALSE))
  expect_silent(minimexico(mexico_abb, rep("black", 32), border_colors = rep("white", 32), 
                        estados_name_colors = rep("yellow", 32), estados_name_cex = .5, 
                        font = "serif"))
})

test_that("minimexico() throws appropriate errors", {
  expect_error(minimexico(c(mexico_abb[-1], "DIF"), 1:32))
  expect_error(minimexico(mexico_abb[1:10], 1:32))
  expect_error(minimexico(mexico_abb, 1:10))
  expect_error(minimexico(mexico_abb, 1:32, border_colors = 1:10))
  expect_error(minimexico(mexico_abb, 1:32, estados_names = NA))
  expect_error(minimexico(mexico_abb, 1:32, estados_name_colors = 1:10))
})