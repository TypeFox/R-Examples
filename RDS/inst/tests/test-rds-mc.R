# TODO: Add comment
# 
# Author: ianfellows
###############################################################################


library(testthat)
library(RDS)

context("rds-mc.R")

test_that("rds-mc",{
	data(faux)
	tc <- count.transitions(faux,"X")
	expect_identical(tc,
					structure(c(39L, 79L, 82L, 188L), .Dim = c(2L, 2L), .Dimnames = structure(list(
											rgrp = c("blue", "red"), grp = c("blue", "red")), .Names = c("rgrp", 
											"grp")))
	)
	
	expect_equal(get.stationary.distribution(prop.table(tc,1)),
			structure(c(0.303913776110387, 0.696086223889613), 
					status = "MLE gives unique stationary distribution.", 
					.Names = c("blue", "red")))

})
