#
# file: tests/testthat/test_rdist90ci_exact.R
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
# ToDo: is this necessary?
library(decisionSupport)

##############################################################################################
# Test rdist90ci_exact(distribution, n, percentiles, quantiles)
##############################################################################################
context("Checking rdist90ci_exact()")

set.seed(100)
n=10000
tolerance=2/sqrt(n)
##############################################################################################
## test distribution="const"
test_that("A constant is generated correctly from the 0.05 and 0.95 quantiles", {
  n=100
  lower=4.5
  upper=lower
  x<-rdist90ci_exact(distribution="const", n=n, lower=lower, upper=upper)
  expect_equal(x, rep(lower,n))
  expect_equal(mean(x), lower)
  expect_equal(sd(x), 0)
})
##############################################################################################
# test distribution="norm"
test_that("Standard normal distribution is generated correctly from the 0.05 and 0.95 quantiles", {
  mean=0
  sd=1
  percentiles=c(0.05, 0.95)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  x<-rdist90ci_exact(distribution="norm", n=n, lower=quantiles[[1]], upper=quantiles[[2]])
  expect_equal(mean(x), mean, tolerance=tolerance)
  expect_equal(sd(x), sd, tolerance=tolerance)
})
test_that("A normal distribution is generated correctly from the 0.05 and 0.95 quantiles (1)", {
  mean=3
  sd=0.5
  percentiles=c(0.05, 0.95)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  x<-rdist90ci_exact(distribution="norm", n=n, lower=quantiles[[1]], upper=quantiles[[2]])
  expect_equal(mean(x), mean, tolerance=tolerance)
  expect_equal(sd(x), sd, tolerance=tolerance)
})
test_that("A normal distribution is generated correctly from the 0.05 and 0.95 quantiles (2)", {
  mean=60000
  sd=24000
  percentiles=c(0.05, 0.95)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  x<-rdist90ci_exact(distribution="norm", n=n, lower=quantiles[[1]], upper=quantiles[[2]])
  expect_equal(mean(x), mean, tolerance=tolerance)
  expect_equal(sd(x), sd, tolerance=tolerance)
})
test_that("A normal distribution is generated correctly from the 0.05 and 0.95 quantiles (3)", {
  lower=20000
  upper=100000
  percentiles=c(0.05, 0.95)
  x<-rdist90ci_exact(distribution="norm", n=n, lower=lower, upper=upper)
  expect_equal(quantile(x,probs=0.05)[["5%"]], lower, tolerance=1.5*tolerance)
  expect_equal(quantile(x,probs=0.95)[["95%"]], upper, tolerance=1.5*tolerance)
})
test_that("A normal distribution cannot be generated from the lower CI value (upper not given)", {
  mean=3
  sd=0.5
  percentiles=c(0.05)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  expect_error(rdist90ci_exact(distribution="norm", n=n, lower=quantiles[[1]]))
})
test_that("A normal distribution cannot be generated from the lower CI value (upper=NULL)", {
  mean=3
  sd=0.5
  percentiles=c(0.05)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  expect_error(rdist90ci_exact(distribution="norm", n=n, lower=quantiles[[1]], upper=NULL))
})
test_that("A normal distribution cannot be generated from the lower CI value (upper=NA)", {
  mean=3
  sd=0.5
  percentiles=c(0.05)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  expect_error(rdist90ci_exact(distribution="norm", n=n, lower=quantiles[[1]], upper=NA))
})
##############################################################################################
## test distribution="lnorm"
test_that("A log normal distribution is generated correctly from the 0.05 and 0.95 quantiles (1)", {
  meanlog=3
  sdlog=0.5
  percentiles=c(0.05, 0.95)
  quantiles=qlnorm(p=percentiles, meanlog=meanlog, sdlog=sdlog)
  x<-rdist90ci_exact(distribution="lnorm", n=n, lower=quantiles[[1]], upper=quantiles[[2]])
  expect_equal(mean(log(x)), meanlog, tolerance=tolerance)
  expect_equal(sd(log(x)), sdlog, tolerance=tolerance)
})
test_that("A log normal distribution is generated correctly from the 0.05 and 0.95 quantiles (2)", {
  lower=3
  upper=6
  percentiles=c(0.05, 0.95)
  x<-rdist90ci_exact(distribution="lnorm", n=n, lower=lower, upper=upper)
  expect_equal(quantile(x,probs=0.05)[["5%"]], lower, tolerance=1.5*tolerance)
  expect_equal(quantile(x,probs=0.95)[["95%"]], upper, tolerance=1.5*tolerance)
})
test_that("A log normal distribution is generated correctly from the 0.05 and 0.95 quantiles (3)", {
  lower=20000
  upper=100000
  percentiles=c(0.05, 0.95)
  x<-rdist90ci_exact(distribution="lnorm", n=n, lower=lower, upper=upper)
  expect_equal(quantile(x,probs=0.05)[["5%"]], lower, tolerance=1.5*tolerance)
  expect_equal(quantile(x,probs=0.95)[["95%"]], upper, tolerance=1.5*tolerance)
})
##############################################################################################
## test distribution="unif"
test_that("A uniform distribution is generated correctly from the 0.05 and 0.95 quantiles (1)", {
  min=-3
  max=10
  percentiles=c(0.05, 0.95)
  quantiles=qunif(p=percentiles, min=min, max=max)
  x<-rdist90ci_exact(distribution="unif", n=n, lower=quantiles[[1]], upper=quantiles[[2]])
  expect_equal(mean(c(min(x),max(x))), mean(c(min,max)), tolerance=tolerance)
})
test_that("A uniform distribution is generated correctly from the 0.05 and 0.95 quantiles (2)", {
  lower=3
  upper=6
  percentiles=c(0.05, 0.95)
  x<-rdist90ci_exact(distribution="unif", n=n, lower=lower, upper=upper)
  expect_equal(quantile(x,probs=0.05)[["5%"]], lower, tolerance=1.5*tolerance)
  expect_equal(quantile(x,probs=0.95)[["95%"]], upper, tolerance=1.5*tolerance)
})
test_that("A uniform distribution is generated correctly from the 0.05 and 0.95 quantiles (3)", {
  lower=20000
  upper=100000
  percentiles=c(0.05, 0.95)
  x<-rdist90ci_exact(distribution="unif", n=n, lower=lower, upper=upper)
  expect_equal(quantile(x,probs=0.05)[["5%"]], lower, tolerance=2*tolerance)
  expect_equal(quantile(x,probs=0.95)[["95%"]], upper, tolerance=2*tolerance)
})
