# create data
X = matrix(c(0,0,1,
             0,1,1,
             1,0,1,
             1,1,1), nrow=4, byrow=TRUE)

y = matrix(c(0,
             1,
             1,
             0),
           nrow=4)

# run full function
learn_bp(X, y)

# define context
context("output format")

test_that("output dimensions match", {
  expect_equal( dim(syn0), c(3, 4) )
  expect_equal( dim(syn1), c(4, 1) )
})
