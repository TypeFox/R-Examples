#
# file: tests/testthat/test_random.estimate.R
#
# This file is part of the R-package decisionSupport
# 
# Authors: 
#   Lutz Göhring <lutz.goehring@gmx.de>
#   Eike Luedeling (ICRAF) <eike@eikeluedeling.com>
#
# Copyright (C) 2015 World Agroforestry Centre (ICRAF) 
#	http://www.worldagroforestry.org
# 
# The R-package decisionSupport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# The R-package decisionSupport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with the R-package decisionSupport.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################################
context("Testing random.estimate()")

set.seed(100)
# Number of random numbers to be generated:
n= 10000
tolerance=3/sqrt(n)

test_that("3d - standard normal distribution (correlated) is generated correctly from the 0.05 and 0.95 quantiles", {
	# Read the estimate for revenue and costs:
	profitEstimate<-estimate_read_csv("profit-4.csv")
	# Calculate means from 95%-confidence intervalls:
	meanProductPrice <- mean(c(profitEstimate$marginal["productprice","lower"], profitEstimate$marginal["productprice","upper"]) ) 
	meanCostPrice <- mean(c(profitEstimate$marginal["costprice","lower"], profitEstimate$marginal["costprice","upper"]) ) 
	meanSales <- mean(c(profitEstimate$marginal["sales","lower"], profitEstimate$marginal["sales","upper"]) )
	mean <- c(productprice=meanProductPrice, costprice=meanCostPrice, sales=meanSales)
	# Calculate standard deviations from 95%-confidence intervalls:
	sdProductPrice <- 0.5 * (profitEstimate$marginal["productprice","upper"] - profitEstimate$marginal["productprice","lower"]) / qnorm(0.95)
	sdCostPrice <- 0.5 * (profitEstimate$marginal["costprice","upper"] - profitEstimate$marginal["costprice","lower"]) / qnorm(0.95)
	sdSales <- 0.5 * (profitEstimate$marginal["sales","upper"] - profitEstimate$marginal["sales","lower"]) / qnorm(0.95)
	sd <- c(productprice=sdProductPrice, costprice=sdCostPrice, sales=sdSales)
	# Extract the full correlation matrix:
	cor<-corMat(profitEstimate)
	# Generate the random numbers:
	x<-random(rho=profitEstimate,n=n)
	# Test:
	expect_equal(colMeans(x), mean, tolerance=tolerance)	
	expect_equal(apply(X=x, MARGIN=2, sd), sd, tolerance=tolerance)
	expect_equal(cor(x),cor, tolerance=0.05) 
})
