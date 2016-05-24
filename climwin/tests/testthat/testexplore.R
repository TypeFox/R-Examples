# Test the outcomes of the explore function #

# Test that errors are thrown correctly #
test_that("explore produces parameter errors", {
  expect_error(explore(shape=-1, scale=0.2, loc=0, weightfunc = "W"))
  expect_error(explore(shape=3, scale=-1, loc=0, weightfunc = "W"))
  expect_error(explore(shape=3, scale=0.2, loc=1, weightfunc = "W"))
  expect_error(explore(shape=3, scale=-1, loc=0, weightfunc = "G"))
})