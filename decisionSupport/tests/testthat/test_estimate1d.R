#
# file: tests/testthat/test_estimate1d.R
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
# method: random.estimate1d()
##############################################################################################
context("Testing random.estimate1d()")

set.seed(100)
# Number of random numbers to be generated:
n= 1000
tolerance=3/sqrt(n)

test_that("method=\"calculate\": Warning is generated for deviation from median if supplied", {
  method="calculate"
  lower<-5
  upper<-20
  c_0.95=qnorm(0.95)
  meanlog<-mean( c(log(lower),log(upper)) )
  sdlog<-( meanlog - log(lower) )/c_0.95
  median_exact<-qlnorm(p=0.50,meanlog = meanlog, sdlog = sdlog)
  x<-random(estimate1d("lnorm", lower,upper),n=n, method=method, relativeTolerance=tolerance)  
  x<-random(estimate1d("lnorm", lower,upper, median=NULL),n=n, method=method, relativeTolerance=tolerance)
  expect_warning(x<-random(estimate1d("lnorm", lower,upper, median="mean"),  n=n, method=method, relativeTolerance=tolerance))
  expect_warning(x<-random(estimate1d("lnorm", lower,upper, median=upper-1), n=n, method=method, relativeTolerance=tolerance))
  x<-random(estimate1d("lnorm", lower,upper, median=median_exact), n=n, method=method, relativeTolerance=tolerance)
})
test_that("method=\"fit\": Warning is generated if no median is supplied or for deviation from median (if supplied)", {
  method="fit"
  lower<-5
  upper<-20
  c_0.95=qnorm(0.95)
  meanlog<-mean( c(log(lower),log(upper)) )
  sdlog<-( meanlog - log(lower) )/c_0.95
  median_exact<-qlnorm(p=0.50,meanlog = meanlog, sdlog = sdlog)
  expect_warning(x<-random(estimate1d("lnorm", lower,upper),n=n, method=method, relativeTolerance=tolerance)) 
  expect_warning(x<-random(estimate1d("lnorm", lower,upper, median=NULL),n=n, method=method, relativeTolerance=tolerance))
  expect_warning(x<-random(estimate1d("lnorm", lower,upper, median="mean"),  n=n, method=method, relativeTolerance=tolerance))
  expect_warning(x<-random(estimate1d("lnorm", lower,upper, median=upper-1), n=n, method=method, relativeTolerance=tolerance))
  x<-random(estimate1d("lnorm", lower,upper, median=median_exact), n=n, method=method, relativeTolerance=tolerance)
})
