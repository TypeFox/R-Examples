context("Modularity")

test_that("Modularity returns values (non NA)", {
	set.seed(100);
  test <- matrix(rbinom(100,1,0.5), ncol=10)
  mod.test <- Modularity(test, sims=10, nstarts=10)
  expect_equal(any(is.na(mod.test)), FALSE)
})
