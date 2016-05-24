

library(testthat)
library(RDS)

context("differential-activity.R")

test_that("differential-activity",{
			data(faux)
			expect_equal(structure(c(1.03342283930911, 1), class = "differential.activity.estimate"),
					differential.activity.estimates(faux,"X",weight.type="RDS-II"))
		})