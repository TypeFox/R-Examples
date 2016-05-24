context("Metacommunity")

test_that("Metacommunity is equal to sum of parts", {
	set.seed(100);
	test <- matrix(rbinom(100,1,0.5), ncol=10)
	mc.test <- Metacommunity(test, sims=100, method='swap')
	bc.test <- BoundaryClump(test)
	expect_equal(mc.test$Boundary['index'], bc.test['index'])
})
