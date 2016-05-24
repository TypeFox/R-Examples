# TODO: Add comment
# 
# Author: ianfellows
###############################################################################



library(testthat)
library(RDS)

context("rds-II.R")

test_that("rds-II",{
			data(faux)
			
			est <- RDS.II.estimates(rds.data=faux,outcome.variable='X')
			expect_equal(est$estimate,structure(c(0.310334112778454, 0.689665887221546), .Names = c("blue", 
									"red")))
		})
