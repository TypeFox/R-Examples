#
# file: tests/testthat/test_rdistq_fit.R
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
# Test rdistq_fit(distribution, n, percentiles, quantiles)
##############################################################################################
context("Checking rdistq_fit()")
set.seed(100)
n=10000
tolerance=2/sqrt(n)
##############################################################################################
## test distribution="norm"
test_that("Standard normal distribution is generated correctly from the 0.05 and 0.95 quantiles", {
  mean=0
  sd=1
  percentiles=c(0.05, 0.95)
  (quantiles=qnorm(p=percentiles, mean=mean, sd=sd))
  x<-rdistq_fit(distribution="norm", n=n, percentiles=percentiles, quantiles=quantiles)
  expect_equal(mean(x), mean, tolerance=tolerance)
  expect_equal(sd(x), sd, tolerance=tolerance)
})
test_that("A normal distribution is generated correctly from the 0.05 and 0.95 quantiles (1)", {
  mean=3
  sd=0.5
  percentiles=c(0.05, 0.95)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  x<-rdistq_fit(distribution="norm", n=n, percentiles=percentiles, quantiles=quantiles)
  expect_equal(mean(x), mean, tolerance=tolerance)
  expect_equal(sd(x), sd, tolerance=tolerance)
})
test_that("A normal distribution is generated correctly from the 0.05 and 0.95 quantiles, if no warning occurs (1)", {
  mean=60000
  sd=24000
  percentiles=c(0.05, 0.95)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  if( inherits(try(expr=expect_warning(x<-rdistq_fit(distribution="norm", n=n, percentiles=percentiles, quantiles=quantiles))),
               "try-error")){
    expect_equal(mean(x), mean, tolerance=tolerance)
    expect_equal(sd(x), sd, tolerance=tolerance)
  }
})
test_that("A normal distribution is generated correctly from the 0.05, 0.10, 0.17, 0.45 and 0.95 quantiles", {
  mean=3
  sd=0.5
  percentiles=c(0.05, 0.10, 0.17, 0.45, 0.95)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  x<-rdistq_fit(distribution="norm", n=n, percentiles=percentiles, quantiles=quantiles)
  expect_equal(mean(x), mean, tolerance=tolerance)
  expect_equal(sd(x), sd, tolerance=tolerance)
})
test_that("A normal distribution cannot be generated from only one quantile, here 0.05", {
  mean=3
  sd=0.5
  percentiles=c(0.05)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  expect_error(rdistq_fit(distribution="norm", n=n, percentiles=percentiles, quantiles=quantiles)
  )
})
##############################################################################################
## test distribution="lnorm"
test_that("A log normal distribution is generated correctly from the 0.05 and 0.95 quantiles (1)", {
  meanlog=3
  sdlog=0.5
  percentiles=c(0.05, 0.95)
  quantiles=qlnorm(p=percentiles, meanlog=meanlog, sdlog=sdlog)
  x<-rdistq_fit(distribution="lnorm", n=n, percentiles=percentiles, quantiles=quantiles)
  expect_equal(mean(log(x)), meanlog, tolerance=tolerance)
  expect_equal(sd(log(x)), sdlog, tolerance=tolerance)
})
test_that("A log normal distribution is generated correctly from the 0.05 and 0.95 quantiles (2)", {
  lower=3
  upper=6
  quantiles<-c("5%"=lower,"95%"=upper)
  percentiles=c(0.05, 0.95)
  x<-rdistq_fit(distribution="lnorm", n=n, percentiles=percentiles, quantiles=quantiles)
  expect_equal(quantile(x,probs=percentiles), quantiles, tolerance=1.5*tolerance)
})
test_that("A log normal distribution is generated correctly from the 0.05 and 0.95 quantiles (3)", {
  lower=20000
  upper=100000
  quantiles<-c("5%"=lower,"95%"=upper)
  percentiles=c(0.05, 0.95)
  x<-rdistq_fit(distribution="lnorm", n=n, percentiles=percentiles, quantiles=quantiles)
  expect_equal(quantile(x,probs=percentiles), quantiles, tolerance=1.5*tolerance)
})
##############################################################################################
## test distribution="unif"
test_that("A uniform distribution is generated correctly from the 0.05 and 0.95 quantiles (1)", {
  min=-3
  max=10
  percentiles=c(0.05, 0.95)
  quantiles=qunif(p=percentiles, min=min, max=max)
  x<-rdistq_fit(distribution="unif", n=n, percentiles=percentiles, quantiles=quantiles)
  expect_equal(mean(c(min(x),max(x))), mean(c(min,max)), tolerance=tolerance)
})
test_that("A uniform distribution is generated correctly from the 0.05 and 0.95 quantiles (2)", {
  lower=3
  upper=6
  quantiles<-c("5%"=lower,"95%"=upper)
  percentiles=c(0.05, 0.95)
  x<-rdistq_fit(distribution="unif", n=n, percentiles=percentiles, quantiles=quantiles)
  expect_equal(quantile(x,probs=percentiles), quantiles, tolerance=1.5*tolerance)
})
test_that("A uniform distribution is generated correctly from the 0.05 and 0.95 quantiles (3)", {
  lower=20000
  upper=100000
  quantiles<-c("5%"=lower,"95%"=upper)
  percentiles=c(0.05, 0.95)
  x<-rdistq_fit(distribution="unif", n=n, percentiles=percentiles, quantiles=quantiles)
  expect_equal(quantile(x,probs=percentiles), quantiles, tolerance=1.5*tolerance)
})

