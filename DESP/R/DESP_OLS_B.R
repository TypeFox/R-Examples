# DESP/R/DESP_OLS.R by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License (version 3) as published by
#  the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

DESP_OLS_B <-
function(X,SPC) {
  # estimation of B by ordinary least squares
  # the observations of the data matrix X are assumed to have zero mean
  # we compute an estimator of each line for every selected (non-zero) coefficient by OLS

  .Call("selective_OLS_B", X, SPC)
}
