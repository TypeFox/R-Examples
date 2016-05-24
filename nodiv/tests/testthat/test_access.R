context("distrib_data access")

data(coquettes)

test_that("occurrences", {
  expect_equal(sum(occurrences(coquettes, 12)), 1765)
  expect_equal(length(occurrences(coquettes, 12)), 28)
  expect_equal(length(occurrences(coquettes, c(12,21))), 2)
  expect_equal(class(occurrences(coquettes, c("Heliangelus_strophianus", "Discosura_popelairii"))), "list")
  expect_equal(as.numeric(sapply(occurrences(coquettes, c("Heliangelus_strophianus", "Discosura_popelairii")), length)), c(5,2))
  expect_named(occurrences(coquettes, c("Heliangelus_strophianus", "Discosura_popelairii")))
  expect_equal(class(occurrences(coquettes, c(13,19), "logical")), "matrix")
  expect_equal(sum(occurrences(coquettes, c(13,19), "logical")), 47)
  expect_equal(dim(occurrences(coquettes, c(13,19), "logical")), c(154, 2))  
  expect_equal(typeof(occurrences(coquettes, c(13,19), "logical")), "logical")
  expect_equal(head(occurrences(coquettes, c("Heliangelus_strophianus", "Discosura_popelairii"), "names")[[1]]), c("186", "239", "256", "335", "349"))
})

test_that("assemblage", {
  
})