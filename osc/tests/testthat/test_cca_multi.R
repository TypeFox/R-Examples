context("CCA Multi")

test_that("cca works for population", {
  data(population)
  expect_equivalent(length(table(cca(population, s=1, mode=1))), 429)
})

test_that("cca works for population mode 3", {
  data(population)
  expect_equivalent(length(table(cca(population, s=10, mode=3))), 15)
})

