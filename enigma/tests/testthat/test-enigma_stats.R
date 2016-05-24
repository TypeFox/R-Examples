# tests for enigma_stats fxn in taxize
context("enigma_stats")

cbase <- 'com.crunchbase.info.companies.acquisition'

test_that("enigma_stats works correctly with varchar column", {
  skip_on_cran()
  
  res1 <- enigma_stats(dataset=cbase, select='acquired_month')
  expect_is(res1, "enigma_stats")
  expect_true(res1$success)
  expect_is(res1$datapath, "character")
  expect_is(res1$info, "list")
  expect_is(res1$result[[1]], "data.frame")
})

test_that("enigma_stats works correctly with numeric column", {
  skip_on_cran()
  
  res2 <- enigma_stats(dataset=cbase, select='price_amount')
  expect_is(res2, "enigma_stats")
  expect_true(res2$success)
  expect_is(res2$datapath, "character")
  expect_is(res2$info, "list")
  expect_is(res2$result[[1]], "data.frame")
})

test_that("enigma_stats works correctly with date column", {
  skip_on_cran()
  
  pakistan <- 'gov.pk.secp.business-registry.all-entities'
  res3 <- enigma_stats(dataset=pakistan, select='registration_date')
  expect_is(res3, "enigma_stats")
  expect_true(res3$success)
  expect_is(res3$datapath, "character")
  expect_is(res3$info, "list")
  expect_is(res3$result[[1]], "data.frame")
})
