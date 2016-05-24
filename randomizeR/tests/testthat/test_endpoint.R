#######################################################################
# --------------------------------------------------------------------#
# Tests for objects from the endpoint class and associated functions  #
# --------------------------------------------------------------------#
#######################################################################
context("Endpoints")

# Normally distributed endpoints
test_that("normEnd returns valid object", {
	mu <-rnorm(2); sigma <- sample(1:10,2, replace=T)
    expect_is(normEndp(mu,sigma), "endpoint")
    
	sigma <- -sample(1:10,2, replace=T)
	expect_error(normEndp(mu,sigma))
	
	mu <-rnorm(3); sigma <- sample(1:10,2, replace=T)
	expect_error(normEndp(mu,sigma))
	
	mu <-rnorm(2); sigma <- sample(1:10,4, replace=T)
	expect_error(normEndp(mu,sigma))
  }
)