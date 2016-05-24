

library(testthat)
library(RDS)

context("homophily.R")

test_that("homophily",{
			data(fauxmadrona)
			res <- homophily.estimates(fauxmadrona,outcome.variable="disease",
					N=1000,weight.type="RDS-II")
			expect_equal(structure(1.60397806961622, class = "table"),res@estimate)
		})


