# test_e2r.R

# load the decompr package
library(rddtools)

# load the example data set
data(house)

# create rdd_data sets
rd<- rdd_data(x=house$x, y=house$y, cutpoint=0)
rd2 <- rdd_data(x=x, y=y, data=house, cutpoint=0)

# define context
context("input ambivalence")

test_that("is rd equal to rd2?", {
  expect_equal( rd, rd2)
}
)

# define context
context("rd: output format")

test_that("rd: output dimensions match", {
  expect_equal( dim(rd), c(6558, 2) )
})

test_that("rd: output values match", {
  expect_equal( rd[1   ,1],  0.1049 )
  expect_equal( rd[1   ,2],  0.581  )
  expect_equal( rd[4   ,1],  0.0868 )
  expect_equal( rd[4   ,2],  0.5846 )
  expect_equal( rd[6558,1], -0.1982 )
  expect_equal( rd[6558,2],  0.802  )
})
