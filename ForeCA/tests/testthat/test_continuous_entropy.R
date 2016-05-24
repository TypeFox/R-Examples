context("continuous_entropy")

# U(-1, 1)
my_density <- function(x){
  dunif(x, -1, 1)
}

test_that("continuous_entropy computes the right entropy for uniform distribution", {
  # entropy of uniform equals log(b - a)
  for (aa in c(1, 3, 5)) {
    expect_equal(log2(aa - 0), 
                 c(continuous_entropy(function(x) dunif(x, 0, aa), 
                                      lower = 0, upper = aa)))
  }
})

test_that("throws error if integration limits are not specified",{
  expect_error(continuous_entropy(my_density))
})

test_that("throws error if density does not integrate to 1",{
  expect_error(continuous_entropy(my_density, -0.5, 1))
})

test_that("throws error if density function is negative",{
  expect_error(continuous_entropy(function(x) x, -0.5, 1))
})
