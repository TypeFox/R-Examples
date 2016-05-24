# test_nrca.R

# load the decompr packages
library(decompr)
library(gvc)

# load the example data set
data(leather)

# define context
context("leontief long format")

# create a leontief long decomposed data set
l <- decomp(inter,
            final,
            countries,
            industries,
            out)

# apply the nrca to leontief
lnrca <- nrca( l )

test_that("output size matches", {
  expect_equal( length(lnrca), 9 )
})

test_that("output values match", {
  expect_equal( lnrca[1], 1.2676927, tolerance = .002)
  expect_equal( lnrca[9], 1.9843574, tolerance = .002)
})

# define short context
context("leontief matrix format")

do <- load_tables_vectors(inter,
                          final,
                          countries,
                          industries,
                          out)

lm <- leontief(do, long=FALSE)

lm_nrca <- nrca(lm)

test_that("output size matches", {
  expect_equal( length(lm_nrca), 9 )
})

test_that("output values match", {
  expect_equal( lm_nrca[1], 1.267693, tolerance = .02, check.attributes=FALSE)
  expect_equal( lm_nrca[9], 1.984357, tolerance = .02, check.attributes=FALSE)
})

