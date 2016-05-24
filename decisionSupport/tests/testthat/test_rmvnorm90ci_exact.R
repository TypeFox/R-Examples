#
# file: tests/testthat/test_rmvnorm90ci_exact.R
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
# Test rmvnorm90ci_exact(distribution, n, percentiles, quantiles)
##############################################################################################
context("Checking rmvnorm90ci_exact()")

set.seed(100)
n=10000
tolerance=2/sqrt(n)

test_that("1d - standard normal distribution is generated correctly from the 0.05 and 0.95 quantiles", {
  mean=0
  sd=1
  percentiles=c(0.05, 0.95)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  x<-rmvnorm90ci_exact(n=n, lower=quantiles[[1]], upper=quantiles[[2]], correlationMatrix=1)
  expect_equal(mean(x), mean, tolerance=tolerance)
  expect_equal(sd(x), sd, tolerance=tolerance)
})

test_that("2d - standard normal distribution (uncorrelated) is generated correctly from the 0.05 and 0.95 quantiles", {
  mean=c(0,0)
  sd=c(1,1)
  cor=diag(c(1,1))
  percentiles=c(0.05, 0.95)
  quantiles<-qnorm(p=percentiles, mean=mean[[1]], sd=sd[[1]])
  quantiles=cbind(quantiles, quantiles, deparse.level=0)
  x<-rmvnorm90ci_exact(n=n, lower=quantiles[1,], upper=quantiles[2,], correlationMatrix=cor)
  expect_equal(colMeans(x), mean, tolerance=tolerance)
  expect_equal(apply(X=x, MARGIN=2, sd), sd, tolerance=tolerance)
  expect_equal(cor(x),cor, tolerance=0.05)
})

test_that("3d - standard normal distribution (uncorrelated) is generated correctly from the 0.05 and 0.95 quantiles", {
  mean=c(0,0,0)
  sd=c(1,1,1)
  cor=diag(c(1,1,1))
  percentiles=c(0.05, 0.95)
  quantiles<-qnorm(p=percentiles, mean=mean[[1]], sd=sd[[1]])
  quantiles=cbind(quantiles, quantiles, quantiles, deparse.level=0)
  x<-rmvnorm90ci_exact(n=n, lower=quantiles[1,], upper=quantiles[2,], correlationMatrix=cor)
  expect_equal(colMeans(x), mean, tolerance=tolerance)
  expect_equal(apply(X=x, MARGIN=2, sd), sd, tolerance=tolerance)
  expect_equal(cor(x),cor, tolerance=0.05)
})

test_that("3d - standard normal distribution (correlated) is generated correctly from the 0.05 and 0.95 quantiles", {
  mean=c(0,0,0)
  sd=c(1,1,1)
  cor<-         t(c( 1,     0.5,    -0.5))
  cor<-rbind(cor, c( 0.5,   1,       0))
  cor<-rbind(cor, c(-0.5,   0,       1))
  percentiles=c(0.05, 0.95)
  quantiles<-qnorm(p=percentiles, mean=mean[[1]], sd=sd[[1]])
  quantiles=cbind(quantiles, quantiles, quantiles, deparse.level=0)
  x<-rmvnorm90ci_exact(n=n, lower=quantiles[1,], upper=quantiles[2,], correlationMatrix=cor)
  expect_equal(colMeans(x), mean, tolerance=tolerance)
  expect_equal(apply(X=x, MARGIN=2, sd), sd, tolerance=tolerance)
  expect_equal(cor(x),cor, tolerance=0.05) 
})

test_that("A 3-d normal distribution (correlated) is generated correctly from the 0.05 and 0.95 quantiles (1)", {
  mean=c(3,-5, 0.5)
  sd=c(0.5, 1.5, 2)
  cor<-         t(c( 1,     0.5,    -0.5))
  cor<-rbind(cor, c( 0.5,   1,       0))
  cor<-rbind(cor, c(-0.5,   0,       1))
  percentiles=c(0.05, 0.95)
  quantiles<-qnorm(p=percentiles, mean=mean[[1]], sd=sd[[1]])
  quantiles<-cbind(quantiles, qnorm(p=percentiles, mean=mean[[2]], sd=sd[[2]]), deparse.level=0)
  quantiles<-cbind(quantiles, qnorm(p=percentiles, mean=mean[[3]], sd=sd[[3]]), deparse.level=0)
  x<-rmvnorm90ci_exact(n=n, lower=quantiles[1,], upper=quantiles[2,], correlationMatrix=cor)
  expect_equal(colMeans(x), mean, tolerance=tolerance)
  expect_equal(apply(X=x, MARGIN=2, sd), sd, tolerance=tolerance)
  expect_equal(cor(x),cor, tolerance=0.05) 
})

test_that("A 3-d normal distribution (correlated) is generated correctly from the 0.05 and 0.95 quantiles (2)", {
  mean=c(60000, 4000, 5)
  sd=c(24000, 7000, 2)
  cor<-         t(c( 1,     0.5,    -0.5))
  cor<-rbind(cor, c( 0.5,   1,       0))
  cor<-rbind(cor, c(-0.5,   0,       1))
  percentiles=c(0.05, 0.95)
  quantiles<-qnorm(p=percentiles, mean=mean[[1]], sd=sd[[1]])
  quantiles<-cbind(quantiles, qnorm(p=percentiles, mean=mean[[2]], sd=sd[[2]]), deparse.level=0)
  quantiles<-cbind(quantiles, qnorm(p=percentiles, mean=mean[[3]], sd=sd[[3]]), deparse.level=0)
  x<-rmvnorm90ci_exact(n=n, lower=quantiles[1,], upper=quantiles[2,], correlationMatrix=cor)
  expect_equal(colMeans(x), mean, tolerance=tolerance)
  expect_equal(apply(X=x, MARGIN=2, sd), sd, tolerance=tolerance)
  expect_equal(cor(x),cor, tolerance=0.05) 
})

test_that("A 3-d normal distribution (correlated) is generated correctly from the 0.05 and 0.95 quantiles (3)", {
  lower=c(20000, -10000, 30000)
  upper=c(100000, 30000, 80000)
  cor<-         t(c( 1,     0.5,    -0.5))
  cor<-rbind(cor, c( 0.5,   1,       0))
  cor<-rbind(cor, c(-0.5,   0,       1))
  percentiles=c(0.05, 0.95)
  x<-rmvnorm90ci_exact(n=n, lower=lower, upper=upper, correlationMatrix=cor)
  expect_equal(apply(X=x,MARGIN=2,FUN=quantile,probs=0.05), lower, tolerance=1.5*tolerance)
  expect_equal(apply(X=x,MARGIN=2,FUN=quantile,probs=0.95), upper, tolerance=1.5*tolerance)
  expect_equal(cor(x),cor, tolerance=0.05) 
})

test_that("A normal distribution cannot be generated from the lower CI value (upper not given)", {
  mean=3
  sd=0.5
  percentiles=c(0.05)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  expect_error(rmvnorm90ci_exact(n=n, lower=quantiles[[1]], correlationMatrix=1))
})

test_that("A normal distribution cannot be generated from the lower CI value (upper=NULL)", {
  mean=3
  sd=0.5
  percentiles=c(0.05)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  expect_error(rmvnorm90ci_exact(n=n, lower=quantiles[[1]], correlationMatrix=1, upper=NULL))
})

test_that("A normal distribution cannot be generated from the lower CI value (upper=NA)", {
  mean=3
  sd=0.5
  percentiles=c(0.05)
  quantiles=qnorm(p=percentiles, mean=mean, sd=sd)
  expect_error(rmvnorm90ci_exact(n=n, lower=quantiles[[1]], correlationMatrix=1, upper=NA))
})

test_that("All diagonal elements of correlationMatrix must be equal to 1", {
	mean=c(60000, 4000, 5)
	sd=c(24000, 7000, 2)
	cor<-         t(c( 1,     0.5,    -0.5))
	cor<-rbind(cor, c( 0.5,   1,       0))
	cor<-rbind(cor, c(-0.5,   0,       0.5))
	percentiles=c(0.05, 0.95)
	quantiles<-qnorm(p=percentiles, mean=mean[[1]], sd=sd[[1]])
	quantiles<-cbind(quantiles, qnorm(p=percentiles, mean=mean[[2]], sd=sd[[2]]), deparse.level=0)
	quantiles<-cbind(quantiles, qnorm(p=percentiles, mean=mean[[3]], sd=sd[[3]]), deparse.level=0)
	expect_error(rmvnorm90ci_exact(n=n, lower=quantiles[1,], upper=quantiles[2,], correlationMatrix=cor))
})
