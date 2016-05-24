library('testthat')
library('lookupTable')

df.input <- cars
response <- 'dist'
feature.boundaries <- list(c(-Inf, 5, 10, 15, 20, 25, Inf))
features.con <- c('speed')
dist.table <- lookupTable(df.input, response, feature.boundaries, features.con)
df.test <- data.frame(speed = c(2, 23, 41, 5, 9, 8))

test_that("initialize returns a lookupTable object",{
  expect_is(dist.table, 'lookupTable')
})

test_that("predict returns right entry values",{
  expect_equal(round(predict(dist.table, df.test)), c(6, 83, 38, 19, 19, 19))
})

