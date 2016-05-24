# test_i2e.R

# load the decompr package
library(decompr)

# load the example data set
data(leather)

# create a leontief decomposed data set
l <- decomp(inter,
            final,
            countries,
            industries,
            out)

# load the package
library(gvc)

# apply the Importing to Export 
li2e <- i2e( l )

# define context
context("output format")

test_that("output size matches", {
  expect_equal( dim(li2e), c(9, 3) )
})

test_that("output order matches", {
  expect_equal( li2e[1,1], factor(c("Argentina", "Turkey", "Germany"))[1] )
  expect_equal( li2e[4,1], factor(c("Argentina", "Turkey", "Germany"))[2])
  expect_equal( li2e[9,1], factor(c("Argentina", "Turkey", "Germany"))[3])
})

test_that("output values match", {
  expect_equal( li2e[1,3], 0.05295042, tolerance = .002)
  expect_equal( li2e[9,3], 0.1720676, tolerance = .002)
})