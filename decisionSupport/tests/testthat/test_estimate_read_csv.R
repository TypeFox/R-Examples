#
# file: tests/testthat/test_estimate_read_csv.R
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
# Test estimate_read_csv(filename)
##############################################################################################
context("Checking estimate_read_csv()")

test_that("Simple marginal file is read correctly (full check)",{
  profit_1_reference<-list(marginal=NULL,correlation_matrix=NULL)
  class(profit_1_reference)<-"estimate"
  profit_1_reference$marginal<-data.frame( row.names   =c("revenue","costs"),
                                       distribution=c("norm", "norm"),
                                       lower       =c(20000,    10000),
                                       median      =c(NA,       NA),
                                       upper       =c(100000,   70000),
                                       stringsAsFactors=FALSE)
  profit_1_estimate<-estimate_read_csv("profit-1.csv")
  expect_equal(profit_1_estimate, profit_1_reference)
})
test_that("Simple marginal and correlation file are read correctly (full check)",{
  profit_2_reference<-list(marginal=NULL,correlation_matrix=NULL)
  class(profit_2_reference)<-"estimate"
  profit_2_reference$marginal<-data.frame( row.names   =c("revenue","costs"),
                                       distribution=c("norm", "norm"),
                                       lower       =c(20000,    10000),
                                       median      =c(NA,       NA),
                                       upper       =c(100000,   70000),
                                       stringsAsFactors=FALSE)
  profit_2_reference$correlation_matrix<-matrix(c(1, 0.5, 
                                                  0.5, 1), byrow=TRUE, nrow=2, ncol=2,
                                                  dimnames=list(c("revenue","costs"),
                                                                c("revenue","costs")))                                                     
  profit_2_estimate<-estimate_read_csv("profit-2.csv")
  expect_equal(profit_2_estimate, profit_2_reference)
})
test_that("Rows without variable name are dropped",{
  profit_3_reference<-list(marginal=NULL,correlation_matrix=NULL)
  class(profit_3_reference)<-"estimate"
  profit_3_reference$marginal<-data.frame( row.names   =c("revenue","costs"),
                                       distribution=c("norm", "norm"),
                                       lower       =c(20000,    10000),
                                       median      =c(NA,       NA),
                                       upper       =c(100000,   70000),
                                       stringsAsFactors=FALSE)                                                   
  profit_3_estimate<-estimate_read_csv("profit-3.csv")
  expect_equal(profit_3_estimate, profit_3_reference)
})
